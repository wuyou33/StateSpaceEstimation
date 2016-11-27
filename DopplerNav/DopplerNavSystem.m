classdef DopplerNavSystem < handle
    % DopplerNavSystem. Describe Doppler navigation system (navigation using measurement of Doppler shift of radial velocity relative to the Sun).
    
    properties (Access = private)
        dmu; % Doppler Measurement Unit
        timeData;
        initArgs;
        earthEphemeris;
        sunEphemeris;
        
        procNoise;
        observNoise;
        inferenceModel;
    end
    
    methods (Access = public)
        function obj = DopplerNavSystem(dmu, timeData, initArgs, earthEphemeris, sunEphemeris)
            narginchk(5, 5);
            
            if ~isa(dmu, 'DopplerMeasurementUnit')
                error('[ DopplerNavSystem::dmu] dmu must be instance of the DopplerMeasurementUnit')
            end
            
            obj.dmu = dmu;
            obj.timeData = timeData;
            obj.initArgs = initArgs;
            obj.earthEphemeris = earthEphemeris;
            obj.sunEphemeris = sunEphemeris;
        end
        
        function stateMatrix = resolve(this, initialState, initialCov, estimatorType, reportProgress)
            narginchk(5, 5);
            
            this.init(estimatorType);
            
            state  = initialState;
            cov    = initialCov;
            simNum = this.timeData.SimulationNumber;
            
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = zeros(this.inferenceModel.stateDimension, simNum);
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            decompCov = this.updateFilterParams(cov, estimatorType);
            particleSet = this.initParticleSet(estimatorType, state, cov, decompCov);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1; % 1*(i == 1)
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                len = endBlock - startBlock + 1;
                if i == 1; stateMatrix(:, 1) = initialState; end
                
                for j = startSample:len
                    sample = j + startBlock - 1;
                    
                    if reportProgress && mod((sample / (simNum - 1))*100, 10) == 0
                        disp(['Completed: ', num2str((sample / (simNum - 1)) * 100),' %' ]);
                    end
                    
                    [state, cov, decompCov, particleSet] = this.estimate(state, cov, decompCov, estimatorType, particleSet, sample, tEpoch);
                    
                    stateMatrix(:, sample) = state;
                end
                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
        end
        
        function init(this, estimatorType)
            narginchk(2, 2);
            args.type  = 'state';
            args.tag   = 'State estimation for Doppler navigation system (using radial velocity to the Sun)';
            args.model = gssmDoppler('init', this.initArgs);
            
            [this.procNoise, this.observNoise, this.inferenceModel] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
        end
        
        function [state, cov, decompCov, particleSet] = estimate(this, state, cov, decompCov, estimator, particleSet, sample, tEpoch)
            narginchk(8, 8);
            
            observ = this.dmu.dopplerShift(sample);
            this.updateModelParams(tEpoch, sample);
            
            [state, cov, decompCov, particleSet] = this.evaluate(state, cov, decompCov, observ, estimator, particleSet);
        end
        
        function particleSet = initParticleSet(this, estimatorType, state, cov, decompCov)
            if strcmp(estimatorType{1}, 'pf')
                numParticles = 2e3;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                particleSet.weights             = ones(1, numParticles) / numParticles;
            elseif strcmp(estimatorType{1}, 'gspf')
                numParticles = 1e3;
                particleSet.particlesNum   = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 35, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'gmsppf')
                numParticles = 1e3;
                particleSet.particlesNum = numParticles;
                initialParticles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
                % fit a N component GMM to initial state distribution
                particleSet.stateGMM = gaussMixtureModelFit(initialParticles, 25, [eps 1e5], 'sqrt', 1e-20);
            elseif strcmp(estimatorType{1}, 'sppf')
                numParticles = 3e2;
                particleSet.particlesNum        = numParticles;
                particleSet.particles           = chol(cov, 'lower')*randn(this.inferenceModel.stateDimension, numParticles) + cvecrep(state, numParticles);
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
            alpha = 1e-3;   % scale factor (UKF parameter) 1e-3
            beta  = 2;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
            kappa = 0;     % 0 is optimal for state dimension = 2 (UKF parameter)
            
            switch estimatorType{1}
                case 'ukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                case 'srukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    decompCov = chol(cov, 'lower');
                case 'cdkf'
                    this.inferenceModel.spkfParams = sqrt(3); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    this.inferenceModel.spkfParams = sqrt(3); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCov = chol(cov, 'lower');
                case {'sckf', 'fdckf'}
                    decompCov = svdDecomposition(cov);
                case 'pf'
                    this.inferenceModel.resampleThreshold = 1;
                    this.inferenceModel.estimateType      = 'mean';
                case 'gspf'
                    this.inferenceModel.estimateType = 'mean';
                case 'sppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
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
                    this.inferenceModel.ghkfParams = 5; % order of gauss-hermite polynomial
                case 'sghqf'
                    this.inferenceModel.sghkfParams = [4 3]; % order of gauss-hermite polynomial & manner
                otherwise
                    % do nothing by default
            end
        end
    end
    
    methods (Access = private)
        function [state, cov, sCov, particleSetEst] = evaluate(this, state, cov, sCov, observ, estimator, particleSet)
            
            particleSetEst = [];
            
            switch estimator{1}
                case 'ukf'
                    [state, cov, this.procNoise, this.observNoise] = ukf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'srukf'
                    [state, sCov, this.procNoise, this.observNoise] = srukf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'cdkf'
                    [state, cov, this.procNoise, this.observNoise] = cdkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'srcdkf'
                    [state, sCov, this.procNoise, this.observNoise] = srcdkf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'ckf'
                    [state, cov, this.procNoise, this.observNoise] = ckf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'sckf'
                    [state, sCov, this.procNoise, this.observNoise] = sckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'fdckf'
                    [state, sCov, this.procNoise, this.observNoise] = fdckf(state, sCov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'pf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = pf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'gspf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gspf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'sppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = sppf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'gmsppf'
                    [state, particleSetEst, this.procNoise, this.observNoise] = gmsppf(particleSet, this.procNoise, this.observNoise, observ, [], [], this.inferenceModel);
                case 'cqkf'
                    [state, cov, this.procNoise, this.observNoise] = cqkf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'ghqf'
                    [state, cov, this.procNoise, this.observNoise] = ghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                case 'sghqf'
                    [state, cov, this.procNoise, this.observNoise] = sghqf(state, cov, this.procNoise, this.observNoise, observ, this.inferenceModel, [], []);
                otherwise
                    error(strcat('not supported filter type: ', estimator{1}));
            end
        end
        
        function updateModelParams(this, tEpoch, sample)
            narginchk(3, 3);
            
            params(1) = tEpoch;
            params(2) = this.timeData.SampleTime;
            params(3) = this.timeData.Time(sample);
            
            model = this.inferenceModel.model.setParams(this.inferenceModel.model, params, this.getEarthEphemeris(sample), this.getSunEphemeris(sample));
            this.inferenceModel.model = model;
        end
        
        function ephemeris = getEarthEphemeris(this, sample)
            narginchk(2, 2);
            
            ephemeris.x = this.earthEphemeris.x(sample);
            ephemeris.y = this.earthEphemeris.y(sample);
            ephemeris.z = this.earthEphemeris.z(sample);
            ephemeris.vx = this.earthEphemeris.vx(sample);
            ephemeris.vy = this.earthEphemeris.vy(sample);
            ephemeris.vz = this.earthEphemeris.vz(sample);
        end
        
        function ephemeris = getSunEphemeris(this, sample)
            narginchk(2, 2);
            
            ephemeris.x = this.sunEphemeris.x(sample);
            ephemeris.y = this.sunEphemeris.y(sample);
            ephemeris.z = this.sunEphemeris.z(sample);
            ephemeris.vx = this.sunEphemeris.vx(sample);
            ephemeris.vy = this.sunEphemeris.vy(sample);
            ephemeris.vz = this.sunEphemeris.vz(sample);
        end
    end
    
end
