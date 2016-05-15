classdef IntegratedInsSns < handle
    % IntegratedInsSns - integrated inertial navigation system (INS) and
    % satellity navigation system (SNS)
    % provides method to solve navigation problem via INS and SNS (GPS)
    
    properties(Access = private)
        initialCov;
        initialState;
        simulationNumber;
        ins;
        sns;
        time;
        gssmInsSnsInitArgs;
        sampleTime;
    end
    
    methods (Access = public)
        
        %% Create instance of the IntegratedInsSns
        % insInitArgs: structure with following fields:
        %   insInitArgs.accBiasMu                    - accelerometer bias mean
        %   insInitArgs.accBiasSigma                 - accelerometer bias RMS
        %   insInitArgs.gyroBiasMu                   - gyro bias mean
        %   insInitArgs.gyroBiasSigma                - gyro bias RMS
        %   insInitArgs.accelerationInBodyFrame      - acceleration in body frame
        %   insInitArgs.angularVelocityInBodyFrame   - angular velocity in body frame
        %   insInitArgs.simulationNumber             - simulation number
        %   insInitArgs.timeMinutes                  - time in minutes;
        %   insInitArgs.visualize                    - 0 - not visualize, 1 - vizualize
        %   insInitArgs.T_till_current_epoch         - T_till_current_epoch;
        %   insInitArgs.sampleTime                   - sample time;
        %   insInitArgs.accNoiseVar                  - accelerometer noise variance matrix;
        %   insInitArgs.gyroNoiseVar                 - gyro noise variance matrix;
        % snsInitArgs: structure with following fields:
        %   snsInitArgs.Trajectory - trajectroy which was determined via SNS
        %   snsInitArgs.Velocity   - velocity which was determined via SNS
        function obj = IntegratedInsSns(insInitArgs, snsInitArgs, initialAcceleration, initialAngularVelocity, initialQuaternion, time, sampleTime)
            initialCov = [1e2*eye(3) zeros(3, 19); ... distance error
                zeros(3, 3) 1e-2*eye(3) zeros(3, 16); ... velocity error
                zeros(4, 6) 1e-10*eye(4) zeros(4, 12); ... quaterni error
                zeros(3, 10) diag((insInitArgs.accBiasSigma.*insInitArgs.accBiasSigma)) zeros(3, 9); ... acceler bias
                zeros(3, 13) diag((insInitArgs.gyroBiasSigma.*insInitArgs.gyroBiasSigma)) zeros(3, 6); ... gyro bias
                zeros(3, 16) 1e-12*eye(3) zeros(3, 3); ... acceler scale factor
                zeros(3, 19) 1e-12*eye(3) ... acceler scale factor
                ];
            
            obj.time = time;
            obj.initialCov = initialCov;
            initialState = initialCondition() + sqrt(initialCov)*randn(22, 1); % sqrt(diag(initialCov))
            obj.initialState = initialState;
            obj.initialState(7:10) = quaternionNormalize(initialState(7:10));
            obj.ins = initInertialNavigationSystem('init', insInitArgs);
            obj.sns = initSatelliteNavigationSystem('init', snsInitArgs);
            obj.simulationNumber = insInitArgs.simulationNumber;
            
            gssmInsSnsInitArgs.initialParams = [insInitArgs.accBiasMu; ...
                insInitArgs.accBiasSigma; ...
                insInitArgs.gyroBiasMu; ...
                insInitArgs.gyroBiasSigma; ...
                initialAcceleration; ...
                initialAngularVelocity; ...
                initialQuaternion; ...
                insInitArgs.sampleTime;
                time(1)];
            
            gssmInsSnsInitArgs.processNoiseMean              = [insInitArgs.accBiasMu; insInitArgs.gyroBiasMu];
            gssmInsSnsInitArgs.processNoiseCovariance        = [diag(insInitArgs.accBiasSigma.*insInitArgs.accBiasSigma) ...
                zeros(3,3); ...
                zeros(3,3) ...
                diag(insInitArgs.gyroBiasSigma.*insInitArgs.gyroBiasSigma)];
            gssmInsSnsInitArgs.observationNoiseMean          = zeros(6, 1);
            gssmInsSnsInitArgs.observationNoiseCovariance    = [5*5*eye(3) zeros(3,3); zeros(3,3) 0.01*0.01*eye(3)];
            obj.gssmInsSnsInitArgs = gssmInsSnsInitArgs;
            obj.sampleTime = sampleTime;
        end
        
        function stateVector = simulate(this, insInitialStateVect, estimatorType, visualize, satellitePhaseStateTrue)
            args.type  = 'state';
            args.tag   = 'State estimation for loosely coupled Ins & Sns integrated system';
            args.model = gssmInsSns('init', this.gssmInsSnsInitArgs);
            
            [processNoise, observationNoise, inferenceDataSet] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
            
            state = this.initialState;
            covState = this.initialCov;
            
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
                    decompCovState = chol(this.initialCov)';
                case 'cdkf'
                    inferenceDataSet.spkfParams = sqrt(70); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    inferenceDataSet.spkfParams = sqrt(25); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCovState = chol(this.initialCov)';
                case 'sckf'
                    decompCovState = svdDecomposition(this.initialCov);
                otherwise
                    % do nothing by default
            end
            
            insMeasurement = insInitialStateVect;
            insMeasurement(7:10) = quaternionNormalize(insMeasurement(7:10));
            
            tmp  = zeros(22, this.simulationNumber);
            inov = zeros(6, this.simulationNumber);
            stateVector = SatellitePhaseSpace(insMeasurement, this.simulationNumber);
            stateVector.AddPhaseState(insMeasurement, 1);
            
            for i = 2:1:this.simulationNumber
                insMeasurement = this.ins.Simulate(insMeasurement, i, this.sampleTime, this.time(i));
                snsMeasurement = this.sns.Simulate(i);
                observation    = insMeasurement(1:6) - snsMeasurement;
                
                modelParams(1:3)    = inferenceDataSet.model.params(1:3);   % accelerationBiasMu
                modelParams(4:6)    = inferenceDataSet.model.params(4:6);   % accelerationBiasSigma
                modelParams(7:9)    = inferenceDataSet.model.params(7:9);   % gyroBiasMu
                modelParams(10:12)  = inferenceDataSet.model.params(10:12); % gyroBiasSigma
                modelParams(13:15)  = this.ins.GetAcceleration(i);
                modelParams(16:18)  = this.ins.GetAngularVelocity(i);
                modelParams(19:22)  = insMeasurement(7:10);                 % quaternion
                modelParams(23)     = inferenceDataSet.model.params(23);    % sampleTime
                modelParams(24)     = this.time(i);
                
                updModel = inferenceDataSet.model.setParams(inferenceDataSet.model, modelParams);
                inferenceDataSet.model = updModel;
                
                switch estimatorType{1}
                    case 'ukf'
                        [state, covState, processNoise, observationNoise, internalParams] = ukf(state, ...
                            covState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    case 'srukf'
                        [state, decompCovState, processNoise, observationNoise, internalParams] = srukf(state, ...
                            decompCovState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    case 'cdkf'
                        [state, covState, processNoise, observationNoise, internalParams] = cdkf(state, ...
                            covState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    case 'srcdkf'
                        [state, decompCovState, processNoise, observationNoise, internalParams] = srcdkf(state, ...
                            decompCovState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    case 'ckf'
                        [state, covState, processNoise, observationNoise, internalParams] = ckf(state, ...
                            covState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    case 'sckf'
                        [state, decompCovState, processNoise, observationNoise, internalParams] = sckf(state, ...
                            decompCovState, ...
                            processNoise, ...
                            observationNoise, ...
                            observation, ...
                            inferenceDataSet);
                    otherwise
                        error('not supported filter type' + estimatorType{1});
                end
                
                stateVector.AddPhaseState(insCorrection(insMeasurement, state(1:10))', i);
                
                tmp(:,i) = state;
                inov(:, i) = internalParams.inov(1:6);
            end
            
            timeMinutes = this.time / 60;
            
            if (visualize)
                SatelliteOrbitVisualization(iterations(k));
                
                figure();
                    plot(timeMinutes', inov(1:3,:));
                    title('innovation');
                    legend('inov_x', 'inov_y', 'inov_z');
                    grid on;
                    hold on;
                
                figure();
                    plot(timeMinutes', inov(4:6,:));
                    title('innovation');
                    legend('inov_v_x', 'inov_v_y', 'inov_v_z');
                    grid on;
                    hold on;
                
                figure();
                    plot(timeMinutes', stateVector.Trajectory - satellitePhaseStateTrue.Trajectory);
                    title('distance error');
                    legend('x axis', 'y axis', 'z axis');
                    grid on;
                
                figure();
                    plot(timeMinutes', stateVector.Velocity - satellitePhaseStateTrue.Velocity);
                    title('velocity error');
                    legend('x axis', 'y axis', 'z axis');
                    grid on;
                
                figure();
                    plot(timeMinutes', stateVector.Rotation - satellitePhaseStateTrue.Rotation);
                    title('quaternion');
                    legend('q_0', 'q_x', 'q_y', 'q_z');
                    grid on;
            end
        end
    end
    
end