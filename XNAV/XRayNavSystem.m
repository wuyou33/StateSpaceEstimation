classdef XRayNavSystem
    % XRayNavSystem Simulate navigation system based on x-ray sources (Pulsars and Quasars)
    
    properties(Access=private)
        earthEphemeris;
        sunEphemeris;
        xRaySources;
        xRayDetector;
    end
    
    methods
        function obj = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, xRayDetector)
            obj.earthEphemeris = earthEphemeris;
            obj.sunEphemeris = sunEphemeris;
            obj.xRaySources = xRaySources;
            obj.xRayDetector = xRayDetector;
        end
    end
    
    methods(Access = public)
        function stateEstimate = Simulate(this, initialState, initialCov, time, estimatorType, tEpoh, sampleTime, visualize)
            if (nargin < 7); error(' not enough input arguments'); end
            if (nargin < 8); visualize = 0; end;
            
            phase = this.xRayDetector.Simulate(time, visualize);
            
            initArgs.initialParams = [...
                tEpoh, ...
                sampleTime, ...
                time(1)...
            ];
            
            initArgs.xRaySources = this.xRaySources;
            initArgs.earthEphemeris = [this.earthEphemeris.x(1); this.earthEphemeris.y(1); this.earthEphemeris.z(1)];
            initArgs.sunEphemeris = [this.sunEphemeris.x(1); this.sunEphemeris.y(1); this.sunEphemeris.z(1)];
            initArgs.observationNoiseMean = zeros(7, 1);
            initArgs.invPeriods = this.getInvPeriods();
            
            observCov = xRayToaCovariance(this.xRaySources, this.xRayDetector.DetectorArea, this.xRayDetector.TimeBucket, this.xRayDetector.BackgroundPhotnRate);
            initArgs.observationNoiseCovariance = diag(observCov);
            
            args.type  = 'state';
            args.tag   = 'State estimation for loosely coupled Ins & Sns integrated system';
            args.model = gssmXNav('init', initArgs);
            [processNoise, observationNoise, inferenceDataSet] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
            
            state    = initialState;
            covState = initialCov;
            
            switch estimatorType{1}
                case 'ukf'
                    alpha = 1e-3; % scale factor (UKF parameter)
                    beta  = 2;    % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 0.75; % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    inferenceDataSet.spkfParams = [alpha beta kappa];
                case 'srukf'
                    alpha = 1e-3; % scale factor (UKF parameter)
                    beta  = 2;    % 2 is a optimal setting for Gaussian priors (SRUKF parameter)
                    kappa = 0.75; % 0 is optimal for state dimension = 2 (SRUKF parameter)
                    
                    inferenceDataSet.spkfParams = [alpha beta kappa];
                    decompCovState = chol(covState)';
                case 'cdkf'
                    inferenceDataSet.spkfParams = sqrt(3); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    inferenceDataSet.spkfParams = sqrt(3); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCovState = chol(covState)';
                case 'sckf'
                    decompCovState = svdDecomposition(covState);
                otherwise
                    % do nothing by default
            end
            
            simulationNumber = length(time);            
            stateEstimate    = zeros(6, simulationNumber);                      
                                                
            for i = 2:1:simulationNumber
                if mod((i / simulationNumber)*100, 5) == 0
                    disp(['Completed: ', num2str((i / simulationNumber) * 100),' %' ]);
                end
                
                observation = phase(i,:)';
                
                modelParams(1)      = tEpoh;
                modelParams(2)      = sampleTime;
                modelParams(3)      = time(i);
                
                earthEphemerisStep  = [this.earthEphemeris.x(i); this.earthEphemeris.y(i); this.earthEphemeris.z(i)];
                sunEphemerisStep  = [this.sunEphemeris.x(i); this.sunEphemeris.y(i); this.sunEphemeris.z(i)];
                
                updModel = inferenceDataSet.model.setParams(inferenceDataSet.model, ...
                    modelParams, ...
                    inferenceDataSet.model.xRaySources, ...
                    earthEphemerisStep, ...
                    sunEphemerisStep, ...
                    inferenceDataSet.model.invPeriods ...
                );
                inferenceDataSet.model = updModel;
                
                switch estimatorType{1}
                    case 'ukf'
                        [state, covState, processNoise, observationNoise] = ukf(state, covState, processNoise, observationNoise, observation, inferenceDataSet);
                    case 'srukf'
                        [state, decompCovState, processNoise, observationNoise] = srukf(state, decompCovState, processNoise, observationNoise, observation, inferenceDataSet);
                    case 'cdkf'
                        [state, covState, processNoise, observationNoise] = cdkf(state, covState, processNoise, observationNoise, observation, inferenceDataSet);
                    case 'srcdkf'
                        [state, decompCovState, processNoise, observationNoise] = srcdkf(state, decompCovState, processNoise, observationNoise, observation, inferenceDataSet);
                    case 'ckf'
                        [state, covState, processNoise, observationNoise] = ckf(state, covState, processNoise, observationNoise, observation, inferenceDataSet);
                    case 'sckf'
                        [state, decompCovState, processNoise, observationNoise] = sckf(state, decompCovState, processNoise, observationNoise, observation, inferenceDataSet);
                    otherwise
                        error('not supported filter type' + estimatorType{1});
                end
                
                stateEstimate(:, i) = state';
            end
        end
    end
    
    methods (Access=private)
        function invPeriods = getInvPeriods(this)
            dimension = length(this.xRaySources);
            
            invPeriods   = zeros(1, dimension);
            for i = 1:dimension
                x = this.xRaySources(i);
                invPeriods(i) = x.TwoPiOnPeriod;
            end
        end
    
    end
end