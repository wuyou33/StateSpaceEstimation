classdef IntegratedInsSns < handle
    % IntegratedInsSns - integrated inertial navigation system (INS) and
    % satellity navigation system (SNS)
    % provides method to solve navigation problem via INS and SNS (GPS)
    
    properties(Access = private)
        ins;
        sns;
        timeData;
        procNoise;
        observNoise;
        inferenceModel;
        initArgs;
    end
    
    methods (Access = public)
        function obj = IntegratedInsSns(ins, sns, timeData, initArgs)
            obj.ins = ins;
            obj.sns = sns;
            obj.timeData = timeData;
            obj.initArgs = initArgs;
        end
        
        function stateMatrix = simulate(this, initalState, initialCov, insInitialState, estimatorType, visualize)
            args.type  = 'state';
            args.tag   = 'State estimation for loosely coupled Ins & Sns integrated system';
            args.model = gssmInsSns('init', this.initArgs);
            
            [this.procNoise, this.observNoise, this.inferenceModel] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
                        
            state  = initalState;
            cov    = initialCov;
            simNum = this.timeData.SimulationNumber;
                        
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = SatellitePhaseSpace(zeros(10, 1), simNum);
            insErrorEst = zeros(22, simNum);
            
            insState = insInitialState;
            filterParams = {};
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            time = this.timeData.Time;
                                   
            decompCov = this.updateFilterParams(cov, estimatorType);
            
            if strcmp(estimatorType{1}, 'pf')
                numParticles = 5e2;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(1*cov, 'lower')*randn(22, numParticles) + cvecrep(state, numParticles);                
                particleSet.particles(7:10, :)  = quaternionNormalize(particleSet.particles(7:10, :));
                particleSet.weights             = (1 / numParticles)*ones(1, numParticles);
            elseif strcmp(estimatorType{1}, 'sppf')
                numParticles = 2e2;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(1*cov, 'lower')*randn(22, numParticles) + cvecrep(state, numParticles);
                particleSet.particles(7:10, :)  = quaternionNormalize(particleSet.particles(7:10, :));
                particleSet.weights             = (1 / numParticles)*ones(1, numParticles);
                particleSet.particlesCov        = repmat(decompCov, [1 1 numParticles]);
                particleSet.processNoise        = this.procNoise;
                particleSet.observationNoise    = this.observNoise;
            else
                particleSet = [];
            end
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1;
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                len = endBlock - startBlock + 1;
                if i == 1
                    stateMatrix.AddPhaseState(insInitialState, 1);
                    insErrorEst(:, 1) = initalState;
                end
                
                for j = startSample:len
                    sample = j + startBlock - 1;
                    
                    if mod((sample / simNum)*100, 5) == 0
                        disp(['Completed: ', num2str((sample / simNum) * 100),' %' ]);
                    end
                                        
                    insState = this.ins.simulate(insState, sample, tEpoch);
                    snsState = this.sns.getState(sample);                    
                    observ   = insState(1:6) - snsState;   
                    
                    this.updateModelParams(state, time(sample), sample, insState(7:10));                 
                    
                    if sample == startSample
                        controlProc = [zeros(3, 1); zeros(3, 1); [1; zeros(3, 1)]];
                    else
                        controlProc = state(1:10);
                    end
                                                                                
                    [state, cov, decompCov, param, particleSet] = this.resolve(state, cov, decompCov, observ, estimatorType, controlProc, particleSet);
                    
                    insErrorEst(:, sample) = state;
                    correctedState = insCorrection(insState, state(1:10));                    
                    stateMatrix.AddPhaseState(correctedState, sample);
                    insState = correctedState;
                    
                    if ~strcmp(estimatorType{1}, {'pf','sppf'})
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
                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
            
            if (visualize)
                this.visualize(filterParams, insErrorEst);
            end
        end
    end
    
    methods (Access = private)
        function [state, cov, sCov, param, particleSetEst] = resolve(this, state, cov, sCov, observ, estimator, control, particleSet)            
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
                case 'sppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = sppf(particleSet, this.procNoise, this.observNoise, observ, control, [], this.inferenceModel);
                otherwise
                    error('not supported filter type' + estimator{1});
            end            
            state(7:10) = quaternionNormalize(state(7:10));
            if ~isempty(particleSetEst)
                particleSetEst.particles(7:10, :) = quaternionNormalize(particleSetEst.particles(7:10, :));
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
        
        function updateModelParams(this, state, time, sample, quaternion)
            modelParams(1:3)    = this.inferenceModel.model.params(1:3);            % accelerationBiasMu
            modelParams(4:6)    = this.inferenceModel.model.params(4:6);            % accelerationBiasSigma
            modelParams(7:9)    = this.inferenceModel.model.params(7:9);            % gyroBiasMu
            modelParams(10:12)  = this.inferenceModel.model.params(10:12);          % gyroBiasSigma
            modelParams(13:15)  = this.getCorrectedAcceleration(sample, state);     % corrected acceleration
            modelParams(16:18)  = this.getCorrectedAngularVelocity(sample, state);  % corrected angular velocity
            modelParams(19:22)  = quaternion;                                       % quaternion
            modelParams(23)     = this.timeData.SampleTime;                         % sampleTime
            modelParams(24)     = time;
            
            updModel = this.inferenceModel.model.setParams(this.inferenceModel.model, modelParams);
            this.inferenceModel.model = updModel;
        end
        
        function decompCov = updateFilterParams(this, cov, estimatorType)
            decompCov = [];
            switch estimatorType{1}
                case 'ukf'
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                case 'srukf'
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    decompCov = chol(cov, 'lower');
                case 'cdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    this.inferenceModel.spkfParams = sqrt(25); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCov = chol(cov, 'lower');
                case {'sckf', 'fdckf'}
                    decompCov = svdDecomposition(cov);
                case 'pf'
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'sppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                otherwise
                    % do nothing by default
            end
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
    end
end