classdef XRayNavSystem < handle
    % XRayNavSystem Simulate navigation system based on X-Ray sources (Pulsars and Quasars)
    
    properties(Access = private)
        earthEphemeris;
        sunEphemeris;
        xRaySources;
        xRayDetector;
        timeData;
        initArgs;
        procNoise;
        observNoise;
        inferenceModel;
    end
    
    properties (Dependent, Access = public)
        SampleTime;
    end
    
    methods
        function obj = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, timeData, initArgs, xRayDetector)
            narginchk(6, 6);
            if ~isa(timeData, 'TimeExt'); error('[ XRayNavSystem ] timeData should be instance of the TimeExt'); end;
            
            obj.earthEphemeris = earthEphemeris;
            obj.sunEphemeris = sunEphemeris;
            obj.xRaySources = xRaySources;
            obj.xRayDetector = xRayDetector;
            obj.timeData = timeData;
            obj.initArgs = initArgs;
        end
    end
    
    methods(Access = public)
        function stateMatrix = resolve(this, initialState, initialCov, estimatorType, reportProgress)
            narginchk(5, 5);
            
            this.init(estimatorType);
            
            state  = initialState;
            cov    = initialCov;
            simNum = this.timeData.SimulationNumber;
            
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = zeros(this.inferenceModel.stateDimension, simNum);
            
            filterParams = {};
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            decompCov = this.updateFilterParams(cov, estimatorType);
            particleSet = this.initParticleSet(estimatorType, state, cov, decompCov);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1;
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                len = endBlock - startBlock + 1*(i == 1);
                if i == 1
                    stateMatrix(:, 1) = initialState;
                end
                
                for j = startSample:len
                    sample = j + startBlock - 1;
                    
                    if reportProgress && mod((sample / (simNum - 1))*100, 10) == 0
                        disp(['Completed: ', num2str((sample / (simNum - 1)) * 100),' %' ]);
                    end
                    
                    [state, cov, decompCov, particleSet, param] = this.estimate(state, cov, decompCov, estimatorType, particleSet, sample, tEpoch);
                    
                    stateMatrix(:, sample) = state;
                    
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
                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
        end
        
        function [state, cov, decompCov, particleSet, param] = estimate(this, state, cov, decompCov, estimator, particleSet, sample, tEpoch)
            narginchk(8, 8);
            
            observ = this.xRayDetector.getTOA(sample);
            this.updateModelParams(tEpoch, sample);
            
            [state, cov, decompCov, particleSet, param] = this.evaluate(state, cov, decompCov, observ, estimator, particleSet);
        end
        
        function init(this, estimatorType)
            narginchk(2, 2);
            
            args.type  = 'state';
            args.tag   = 'State estimation for X-Ray navigation system';
            args.model = gssmXNav('init', this.initArgs);
            [this.procNoise, this.observNoise, this.inferenceModel] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
        end
        
        function particleSet = initParticleSet(this, estimatorType, state, cov, decompCov)
            narginchk(5, 5);
            
            if strcmp(estimatorType{1}, 'pf')
                numParticles = 1e4;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                particleSet.weights             = ones(1, numParticles) / numParticles;
            elseif strcmp(estimatorType{1}, 'gspf')
                numParticles = 2e2;
                particleSet.particlesNum   = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 5, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'gmsppf')
                numParticles = 3e3;
                particleSet.particlesNum = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 9, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'sppf')
                numParticles = 1e2;
                particleSet.particlesNum  = numParticles;
                
                if ~strcmp(this.inferenceModel.model.processNoise.covarianceType, 'sqrt')
                    error(' [initParticleSet::this.inferenceModel.model.processNoise.covarianceType] should have type "sqrt"');
                else
                    particleSet.particles  =  chol(cov) * randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                    % this.inferenceModel.model.processNoise.covariance | chol(this.inferenceModel.model.processNoise.covariance, 'lower')
                end
                
                particleSet.weights             = ones(1, numParticles) / numParticles;
                particleSet.particlesCov        = repmat(decompCov, [1 1 numParticles]);
                particleSet.processNoise        = this.procNoise;
                particleSet.observationNoise    = this.observNoise;
            else
                particleSet = [];
            end
        end
        
        function decompCov = updateFilterParams(this, cov, estimatorType)
            narginchk(3, 3);
            decompCov = [];
            
            alpha = 0.075;   % scale factor (UKF parameter) 1e-3
            beta  = 2.75;    % 2 is a optimal setting for Gaussian priors (UKF parameter)
            kappa = 0.1;     % 0 is optimal for state dimension = 2 (UKF parameter)
            
            switch estimatorType{1}
                case 'ukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                case 'srukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    decompCov = chol(cov, 'lower');
                case 'cdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCov = chol(cov, 'lower');
                case {'sckf', 'fdckf'}
                    decompCov = svdDecomposition(cov); % chol(cov, 'lower');
                case 'pf'
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gspf'
                    this.inferenceModel.estimateType = 'mean';
                case 'sppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 0.1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gmsppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'cqkf'
                    this.inferenceModel.cqkfParams = 9; % order of laguerre polynomial
                case 'ghqf'
                    this.inferenceModel.ghkfParams = 4; % order of gauss-hermite polynomial
                case 'sghqf'
                    this.inferenceModel.sghkfParams = [5 3]; % order of gauss-hermite polynomial & manner
                otherwise
                    % do nothing by default
            end
        end
    end
    
    methods (Access = private)
        function [state, cov, sCov, particleSetEst, param] = evaluate(this, state, cov, sCov, observ, estimator, particleSet)
            param = [];
            particleSetEst = [];
            
            switch estimator{1}
                case 'ukf'
                    [state, cov, this.procNoise, this.observNoise, param] = ukf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'srukf'
                    [state, sCov, this.procNoise, this.observNoise, param] = srukf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'cdkf'
                    [state, cov, this.procNoise, this.observNoise, param] = cdkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'srcdkf'
                    [state, sCov, this.procNoise, this.observNoise, param] = srcdkf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'ckf'
                    [state, cov, this.procNoise, this.observNoise, param] = ckf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'sckf'
                    [state, sCov, this.procNoise, this.observNoise, param] = sckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'fdckf'
                    [state, sCov, this.procNoise, this.observNoise, param] = fdckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'pf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = pf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'gspf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gspf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'sppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = sppf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'gmsppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gmsppf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'cqkf'
                    [state, cov, this.procNoise, this.observNoise, param] = cqkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'ghqf'
                    [state, cov, this.procNoise, this.observNoise, param] = ghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'sghqf'
                    [state, cov, this.procNoise, this.observNoise, param] = sghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'ekf'
                    [state, cov, this.procNoise, this.observNoise, param] = ekf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                otherwise
                    error(strcat('not supported filter type: ', estimator{1}));
            end
        end
        
        function updateModelParams(this, tEpoch, sample)
            modelParams(1) = tEpoch;
            modelParams(2) = this.timeData.SampleTime;
            modelParams(3) = this.timeData.Time(sample);
            
            earthEphemerisStep = [this.earthEphemeris.x(sample); this.earthEphemeris.y(sample); this.earthEphemeris.z(sample)];
            sunEphemerisStep   = [this.sunEphemeris.x(sample); this.sunEphemeris.y(sample); this.sunEphemeris.z(sample)];
            
            updModel = this.inferenceModel.model.setParams(this.inferenceModel.model, ...
                modelParams, ...
                this.inferenceModel.model.xRaySources, ...
                earthEphemerisStep, ...
                sunEphemerisStep, ...
                this.inferenceModel.model.invPeriods, ...
                this.inferenceModel.model.mass, ...
                this.inferenceModel.model.gravityModel, ...
                this.inferenceModel.model.startTime);
            this.inferenceModel.model = updModel;
        end
        
        function invPeriods = getInvPeriods(this)
            invPeriods = getInvPeriods(this.xRaySources);
        end
    end
    
    methods
        function res = get.SampleTime(this)
            res = this.timeData.SampleTime;
        end
    end
    
end
