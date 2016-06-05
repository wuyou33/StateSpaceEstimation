classdef IntegratedInsSns < handle
    % IntegratedInsSns - integrated inertial navigation system (INS) and
    % satellity navigation system (SNS)
    % provides method to solve navigation problem via INS and SNS (GPS)
    
    properties(Access = private)
        initialCov;        
        simulationNumber;
        ins;
        sns;
        time;
        gssmInsSnsInitArgs;
        sampleTime;
        gyroScale;
        accScale;
    end
    
    methods (Access = public)               
        function obj = IntegratedInsSns(insInitArgs, snsInitArgs, initialCov, initialAcceleration, initialAngularVelocity, initialQuaternion, time, sampleTime)
            %% Create instance of the IntegratedInsSns
            %             
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
            %   insInitArgs.accScale                     - scale factor of accelerometer;
            %   insInitArgs.gyroScale                    - scale factor of gyro ;
            %   
            % snsInitArgs: structure with following fields:
            %   snsInitArgs.Trajectory - trajectroy which was determined via SNS
            %   snsInitArgs.Velocity   - velocity which was determined via SNS
            %%
                        
            obj.time = time;
            obj.initialCov = initialCov;                        
            
            obj.ins = initInertialNavigationSystem('init', insInitArgs);
            obj.sns = initSatelliteNavigationSystem('init', snsInitArgs);
            obj.simulationNumber = insInitArgs.simulationNumber;
            obj.gyroScale = insInitArgs.gyroScale;
            obj.accScale  = insInitArgs.accScale;
            
            gssmInsSnsInitArgs.initialParams = [insInitArgs.accBiasMu; ...
                insInitArgs.accBiasSigma; ...
                insInitArgs.gyroBiasMu; ...
                insInitArgs.gyroBiasSigma; ...
                initialAcceleration; ...
                initialAngularVelocity; ...
                initialQuaternion; ...
                sampleTime;
                time(1)];
            
            gssmInsSnsInitArgs.processNoiseMean              = [insInitArgs.accBiasMu; insInitArgs.gyroBiasMu];
            gssmInsSnsInitArgs.processNoiseCovariance        = [diag(insInitArgs.accBiasSigma.*insInitArgs.accBiasSigma) zeros(3,3); ...
                                                                zeros(3,3) diag(insInitArgs.gyroBiasSigma.*insInitArgs.gyroBiasSigma)];
            gssmInsSnsInitArgs.observationNoiseMean          = zeros(6, 1);
            gssmInsSnsInitArgs.observationNoiseCovariance    = [(1e-1)^2*eye(3) zeros(3,3); zeros(3,3) (1e-4)*eye(3)];
            obj.gssmInsSnsInitArgs = gssmInsSnsInitArgs;
            obj.sampleTime = sampleTime;
        end
        
        function stateVector = Simulate(this, insInitialStateVect, estimatorType, visualize, satellitePhaseStateTrue)
            args.type  = 'state';
            args.tag   = 'State estimation for loosely coupled Ins & Sns integrated system';
            args.model = gssmInsSns('init', this.gssmInsSnsInitArgs);
            
            [stateNoise, observNoise, inferenceDataSet] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
            
            state = this.initializeState();
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
            tmp(:, 1)  = state;
            inov = zeros(6, this.simulationNumber-1);
            stateVector = SatellitePhaseSpace(insMeasurement, this.simulationNumber);
            stateVector.AddPhaseState(insMeasurement, 1);
            
            for i = 2:1:this.simulationNumber
                if mod((i / this.simulationNumber)*100, 5) == 0
                    disp(['Completed: ', num2str((i / this.simulationNumber) * 100),' %' ]);
                end
                
                insMeasurement = this.ins.Simulate(insMeasurement, i, this.sampleTime, this.time(i));
                snsMeasurement = this.sns.Simulate(i);
                observation    = insMeasurement(1:6) - snsMeasurement;
                
                modelParams(1:3)    = inferenceDataSet.model.params(1:3);         % accelerationBiasMu
                modelParams(4:6)    = inferenceDataSet.model.params(4:6);         % accelerationBiasSigma
                modelParams(7:9)    = inferenceDataSet.model.params(7:9);         % gyroBiasMu
                modelParams(10:12)  = inferenceDataSet.model.params(10:12);       % gyroBiasSigma                
                modelParams(13:15)  = this.GetCorrectedAcceleration(i, state);    % corrected acceleration                 
                modelParams(16:18)  = this.GetCorrectedAngularVelocity(i, state); % corrected angular velocity
                modelParams(19:22)  = insMeasurement(7:10);                       % quaternion
                modelParams(23)     = inferenceDataSet.model.params(23);          % sampleTime
                modelParams(24)     = this.time(i);
                
                updModel = inferenceDataSet.model.setParams(inferenceDataSet.model, modelParams);
                inferenceDataSet.model = updModel;
                
                switch estimatorType{1}
                    case 'ukf'
                        [state, covState, stateNoise, observNoise, internalParams] = ukf(state, covState, stateNoise, observNoise, observation, inferenceDataSet);
                    case 'srukf'
                        [state, decompCovState, stateNoise, observNoise, internalParams] = srukf(state, decompCovState, stateNoise, observNoise, observation, inferenceDataSet);
                    case 'cdkf'
                        [state, covState, stateNoise, observNoise, internalParams] = cdkf(state, covState, stateNoise, observNoise, observation, inferenceDataSet);
                    case 'srcdkf'
                        [state, decompCovState, stateNoise, observNoise, internalParams] = srcdkf(state, decompCovState, stateNoise, observNoise, observation, inferenceDataSet);
                    case 'ckf'
                        [state, covState, stateNoise, observNoise, internalParams] = ckf(state, covState, stateNoise, observNoise, observation, inferenceDataSet);
                    case 'sckf'
                        [state, decompCovState, stateNoise, observNoise, internalParams] = sckf(state, decompCovState, stateNoise, observNoise, observation, inferenceDataSet);
                    otherwise
                        error('not supported filter type' + estimatorType{1});
                end
                
                stateVector.AddPhaseState(insCorrection(insMeasurement, state(1:10))', i);
                
                tmp(:,i) = state;
                inov(:, i-1) = internalParams.inov(1:6);
            end          
            
%             if (visualize)
%                 this.Visualize(stateVector, satellitePhaseStateTrue);
% %                 this.Visualize(stateVector, satellitePhaseStateTrue, tmp);
%             end
        end
    end
        
    methods (Access = private)
        function acceleration = GetCorrectedAcceleration(this, sample, state)
            accelerationEst = this.ins.GetAcceleration(sample)';
            acceleration = (accelerationEst - state(11:13)) ./ (ones(3, 1) + state(17:19));
        end
        
        function angVelocity = GetCorrectedAngularVelocity(this, sample, state)
            angVelocityEst = this.ins.GetAngularVelocity(sample)';
            angVelocity = (angVelocityEst - state(14:16)) ./ (ones(3, 1) + state(20:22));
        end
        
        function Visualize(this, stateVector, satellitePhaseStateTrue, fullStateVector, inov)
            timeMinutes = this.time / 60;
            
            SatelliteOrbitVisualization(stateVector);
            
            if nargin >= 4
                this.Plot2(timeMinutes', fullStateVector(11:13, :), 'acceleration bias', {'x axis', 'y axis', 'z axis'}, 'acceleration bias, km/sec^2');
                this.Plot2(timeMinutes', fullStateVector(14:16, :), 'gyro bias', {'x axis', 'y axis', 'z axis'}, 'gyro bias, rad/sec');
                this.Plot2(timeMinutes', fullStateVector(17:19, :), 'acceleration scale factor', {'x axis', 'y axis', 'z axis'}, 'acceleration scale factor');
                this.Plot2(timeMinutes', fullStateVector(20:22, :), 'gyro scale factor', {'x axis', 'y axis', 'z axis'}, 'gyro scale factor');
            end
            
            if nargin == 5
                this.Plot2(timeMinutes', inov(1:3,:), 'innovation', {'inov_x', 'inov_y', 'inov_z'}, 'inov distance');
                this.Plot2(timeMinutes', inov(4:6,:), 'innovation', {'inov_v_x', 'inov_v_y', 'inov_v_z'}, 'inov velocity');                
            end
            
            this.Plot2Log(timeMinutes', 1e3*(stateVector.Trajectory - satellitePhaseStateTrue.Trajectory), 'distance error', {'x axis', 'y axis', 'z axis'}, 'distance error, m');
            this.Plot2Log(timeMinutes', 1e3*(stateVector.Velocity - satellitePhaseStateTrue.Velocity), 'velocity error', {'x axis', 'y axis', 'z axis'}, 'velocity error, m/sec');
            
            [yawE, pitchE, rollE] = quat2angle(stateVector.Rotation');
            angleEstimated        = [yawE, pitchE, rollE];
            [yaw, pitch, roll]    = quat2angle(satellitePhaseStateTrue.Rotation');
            angleTrue             = [yaw, pitch, roll];
            this.Plot2Log(timeMinutes', abs(angleEstimated - angleTrue), 'angle rotation error', {'yaw', 'pitch', 'roll'}, 'angle error, rad');
        end
        
        function initialState = initializeState(this)
            initialState = sqrt(this.initialCov)*randn(22, 1); % sqrt(diag(initialCov))
            
            initialState(7) = 1;
            initialState(7:10) = quaternionNormalize(initialState(7:10));
            
            initialState(17:19) = this.accScale + initialState(17:19);
            initialState(20:22) = this.gyroScale + initialState(20:22);
        end
        
        function Plot2(~, x, y, titleText, legendText, yLabelText)            
            plot2(x, y, titleText, legendText, yLabelText, 'time, min');
        end
        
        function Plot2Log(~, x, y, titleText, legendText, yLabelText)  
            semilogy2(x, y, titleText, legendText, yLabelText, 'time, min');            
        end
    end
end