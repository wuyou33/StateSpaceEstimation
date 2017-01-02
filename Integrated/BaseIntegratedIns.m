classdef BaseIntegratedIns < handle
    % BaseIntegratedIns. Provide base methods which allow integrate INS with any other navigation system using tougthly coupled alghorithm
    
    properties (Access = protected)
        ins;
        timeData;
        procNoise;
        observNoise;
        inferenceModel;
        initArgs;
        dimension = 10; % dimension of state space in INS navigation issue
        reconciliationTime = NaN; % re-initialization time for filter [ seconds ]
    end
    
    methods (Access = public)
        function stateMatrix = evaluate(this, initalState, initialCov, insInitialState, estimatorType, visualize, reportProgress)
            narginchk(7, 7);
            
            args.type  = 'state';
            args.tag   = this.tag();
            args.model = gssmInsLooselyCoupled('init', this.initArgs);
            
            [this.procNoise, this.observNoise, this.inferenceModel] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
            
            state  = initalState;
            cov    = initialCov;
            simNum = this.timeData.SimulationNumber;
            
            stateMatrix = SatellitePhaseSpace(zeros(this.dimension, 1), this.ins.SimulationNumber);
            insErrorEst = zeros(this.inferenceModel.stateDimension, simNum);
            
            insState = insInitialState;
            filterParams = {};
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            time = this.timeData.Time;
            
            reconciliationIndexes = this.findReconciliationIndexes(simNum);
            
            decompCov = this.updateFilterParams(cov, estimatorType);
            particleSet = this.buildParticles(estimatorType, state, cov, decompCov);
            
            insIndex = 2;
            for i = 1:num
                startBlock = (i-1)*blockSize + 1; % 1*(i == 1)
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                len = endBlock - startBlock + 1;
                if i == 1
                    stateMatrix.AddPhaseState(insInitialState, 1);
                    insErrorEst(:, 1) = initalState;
                end
                
                for j = startSample:len
                    sample = j + startBlock - 1;
                    
                    if any(reconciliationIndexes == sample)
                        state(10:end)   = initalState(10:end);
                        cov             = initialCov;
                        decompCov       = this.updateFilterParams(cov, estimatorType);
                        particleSet     = this.buildParticles(estimatorType, state, cov, decompCov);
                    end
                    
                    if reportProgress && mod((sample / (simNum - 1))*100, 10) == 0
                        disp(['Completed: ', num2str((sample / (simNum - 1)) * 100),' %' ]);
                    end
                    
                    timeStart = time(sample - 1);
                    timeEnd = time(sample);
                    
                    initialSubSysState = insState;
                    this.updateModelParams(state, time(sample), sample, initialSubSysState);
                    
                    insStateMatrix = this.ins.evaluate(initialSubSysState, timeStart, timeEnd);
                    insState = insStateMatrix(:, end);
                    
                    secondarySystemState = this.evaluateSecondaryState(initialSubSysState, timeStart, timeEnd);
                    
                    observ = insState(1:6) - secondarySystemState(1:6);
                    
                    if sample == startSample && i == 1
                        controlProc = [zeros(3, 1); zeros(3, 1); [1; zeros(3, 1)]];
                    else
                        controlProc = state(1:10);
                    end
                    
                    [state, cov, decompCov, param, particleSet] = this.eval(state, cov, decompCov, observ, estimatorType, controlProc, particleSet);
                    
                    insErrorEst(:, sample) = state;
                    correctedState = insCorrection(insState, state(1:10));
                    insState = correctedState;
                    
                    for k = 1:size(insStateMatrix, 2)
                        stateMatrix.AddPhaseState(insStateMatrix(:, k), insIndex);
                        insIndex = insIndex + 1;
                    end
                    
                    if ~strcmp(estimatorType{1}, {'pf','sppf', 'gspf', 'gmsppf'})
                        filterParams.meanPredictedState(sample, :) = param.meanPredictedState;
                        
                        if sum(strcmp(fieldnames(param), 'predictedStateCov')) == 1
                            filterParams.predictedStateCov(sample, :, :) = param.predictedStateCov;
                        end
                        
                        filterParams.predictedObservMean(sample, :) = param.predictedObservMean;
                        filterParams.inov(sample, :, :) = param.inov;
                        
                        if sum(strcmp(fieldnames(param), 'predictedObservCov')) == 1
                            filterParams.predictedObservCov(sample, :, :) = param.predictedObservCov;
                        end
                        
                        filterParams.filterGain(sample, :, :) = param.filterGain;
                    else
                        filterParams = [];
                    end
                end
                
                startSample = 1;
            end
            
            if visualize; this.visualize(filterParams, insErrorEst); end
        end
    end
    
    methods (Abstract, Access = protected)
        updateFilterParams(this, cov, estimatorType);
        
        evaluateSecondaryState(this, initial, timeStart, timeEnd);
        
        tag(this);
    end
    
    methods (Access = private)
        function reconciliationIndexes = findReconciliationIndexes(this, num)
            if isnumeric(this.reconciliationTime) && this.reconciliationTime > 0
                reconciliationSample = floor(this.reconciliationTime / this.timeData.SampleTime);
                reconciliationIndexes = 1 : reconciliationSample : num;
            else
                reconciliationIndexes = [];
            end
        end
        
        function acceleration = getCorrectedAcceleration(this, sample, state)
            a = this.ins.getAcceleration(sample);
            acceleration = (a - state(11:13)) ./ (ones(3, 1) + state(17:19));
        end
        
        function angVelocity = getCorrectedAngularVelocity(this, sample, state)
            w = this.ins.getAngularVelocity(sample);
            angVelocity = (w - state(14:16)) ./ (ones(3, 1) + state(20:22));
        end
        
        function [state, cov, sCov, param, particleSetEst] = eval(this, state, cov, sCov, observ, estimator, control, particleSet)
            param = [];
            particleSetEst = [];
            
            switch estimator{1}
                case 'ukf'
                    [state, cov, this.procNoise, this.observNoise, param] = ukf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'srukf'
                    [state, sCov, this.procNoise, this.observNoise, param] = srukf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'cdkf'
                    [state, cov, this.procNoise, this.observNoise, param] = cdkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'srcdkf'
                    [state, sCov, this.procNoise, this.observNoise, param] = srcdkf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'ckf'
                    [state, cov, this.procNoise, this.observNoise, param] = ckf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'sckf'
                    [state, sCov, this.procNoise, this.observNoise, param] = sckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'fdckf'
                    [state, sCov, this.procNoise, this.observNoise, param] = fdckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'pf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = pf(particleSet, this.procNoise, this.observNoise, observ, control, [], this.inferenceModel);
                case 'gspf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gspf(particleSet, this.procNoise, this.observNoise, observ, control, [], this.inferenceModel);
                    particleSetEst.stateGMM.mean(7:10, :) = quaternionNormalize(particleSetEst.stateGMM.mean(7:10, :));
                case 'sppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = sppf(particleSet, this.procNoise, this.observNoise, observ, control, [], this.inferenceModel);
                case 'gmsppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gmsppf(particleSet, this.procNoise, this.observNoise, observ, control, [], this.inferenceModel);
                    particleSetEst.stateGMM.mean(7:10, :) = quaternionNormalize(particleSetEst.stateGMM.mean(7:10, :));
                case 'cqkf'
                    [state, cov, this.procNoise, this.observNoise, param] = cqkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'ghqf'
                    [state, cov, this.procNoise, this.observNoise, param] = ghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'sghqf'
                    [state, cov, this.procNoise, this.observNoise, param] = sghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                case 'ekf'
                    [state, cov, this.procNoise, this.observNoise, param] = ekf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, control, []);
                otherwise
                    error(strcat('not supported filter type: ', estimator{1}));
            end
            state(7:10) = quaternionNormalize(state(7:10));
            if ~isempty(particleSetEst) && isfield(particleSetEst, 'particles')
                particleSetEst.particles(7:10, :) = quaternionNormalize(particleSetEst.particles(7:10, :));
            end
        end
        
        function updateModelParams(this, state, time, sample, insState)
            modelParams(1:3)    = this.inferenceModel.model.params(1:3);            % accelerationBiasMu
            modelParams(4:6)    = this.inferenceModel.model.params(4:6);            % accelerationBiasSigma
            modelParams(7:9)    = this.inferenceModel.model.params(7:9);            % gyroBiasMu
            modelParams(10:12)  = this.inferenceModel.model.params(10:12);          % gyroBiasSigma
            modelParams(13:15)  = this.getCorrectedAcceleration(sample, state);     % corrected acceleration
            modelParams(16:18)  = this.getCorrectedAngularVelocity(sample, state);  % corrected angular velocity
            modelParams(19:28)  = insState;                                         % attitude, velocity, quaternion of INS
            modelParams(29)     = this.timeData.SampleTime;                         % sampleTime
            modelParams(30)     = time;
            
            updModel = this.inferenceModel.model.setParams(this.inferenceModel.model, modelParams);
            this.inferenceModel.model = updModel;
        end
        
        function visualize(this, filterParams, insErrorEst)
            if ~isempty(insErrorEst)
                figure();
                subplot(3, 1, 1);
                plot2(this.timeData.RelTime, 1e3*insErrorEst(1:3, :), 'estimation trajectory error', {'x', 'y', 'z'}, 'coordinate error, meter');
                subplot(3, 1, 2);
                plot2(this.timeData.RelTime, 1e3*insErrorEst(4:6, :), 'estimation velocity error', {'x', 'y', 'z'}, 'velocity error, meter / sec');
                subplot(3, 1, 3);
                plot2(this.timeData.RelTime, insErrorEst(7:10, :), 'estimation quaternion error', {'q0', 'q1', 'q2', 'q3'}, 'quaternion error,');
                
                figure();
                subplot(2, 2, 1);
                plot2(this.timeData.RelTime, insErrorEst(11:13, :), 'acceleration bias', {'x', 'y', 'z'}, 'acceleration bias, km / sec^2');
                subplot(2, 2, 2);
                plot2(this.timeData.RelTime, insErrorEst(14:16, :), 'gyro bias', {'x', 'y', 'z'}, 'gyro bias, rad / sec');
                subplot(2, 2, 3);
                plot2(this.timeData.RelTime, insErrorEst(17:19, :), 'acceleration scale', {'x', 'y', 'z'}, 'acceleration scale');
                subplot(2, 2, 4);
                plot2(this.timeData.RelTime, insErrorEst(20:22, :), 'gyro scale', {'x', 'y', 'z'}, 'gyro scale');
            end
            
            if ~isempty(filterParams)
                gain = filterParams.filterGain;
                figure();
                subplot(3, 2, 1);
                plot2(this.timeData.RelTime, gain(:, 1, 1), 'filter gain trajectory', {'x'}, 'filter gain');
                subplot(3, 2, 3);
                plot2(this.timeData.RelTime, gain(:, 2, 2), 'filter gain trajectory', {'y'}, 'filter gain');
                subplot(3, 2, 5);
                plot2(this.timeData.RelTime, gain(:, 3, 3), 'filter gain trajectory', {'z'}, 'filter gain');
                subplot(3, 2, 2);
                plot2(this.timeData.RelTime, gain(:, 4, 4), 'filter gain velocity', {'x'}, 'filter gain');
                subplot(3, 2, 4);
                plot2(this.timeData.RelTime, gain(:, 5, 5), 'filter gain velocity', {'y'}, 'filter gain');
                subplot(3, 2, 6);
                plot2(this.timeData.RelTime, gain(:, 6, 6), 'filter gain velocity', {'z'}, 'filter gain');
            end
        end
        
        function particleSet = buildParticles(this, estimatorType, state, cov, decompCov)
            if strcmp(estimatorType{1}, 'pf')
                numParticles = 2e3;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                particleSet.particles(7:10, :)  = quaternionNormalize(particleSet.particles(7:10, :));
                particleSet.weights             = ones(1, numParticles) / numParticles;
            elseif strcmp(estimatorType{1}, 'gspf')
                numParticles = 1e3;
                particleSet.particlesNum   = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                initialParticles(7:10, :)  = quaternionNormalize(initialParticles(7:10, :));
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 10, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'gmsppf')
                numParticles = 5e1;
                particleSet.particlesNum   = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                initialParticles(7:10, :)  = quaternionNormalize(initialParticles(7:10, :));
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 5, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'sppf')
                numParticles = 5e1;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(cov, 'lower')*randn(22, numParticles) + cvecrep(state, numParticles);
                particleSet.particles(7:10, :)  = quaternionNormalize(particleSet.particles(7:10, :));
                particleSet.weights             = ones(1, numParticles) / numParticles;
                particleSet.particlesCov        = repmat(decompCov, [1 1 numParticles]);
                particleSet.processNoise        = this.procNoise;
                particleSet.observationNoise    = this.observNoise;
            else
                particleSet = [];
            end
        end
    end
    
end
