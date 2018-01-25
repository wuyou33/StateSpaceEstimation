close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

set(0, 'defaultfigurecolor', [1 1 1]);
% ghqf not working due to dimension, too number of points :(
% {'ekf', 'ukf', 'cdkf', 'ckf', 'sckf', 'srukf', 'srcdkf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'pf', 'sppf', 'gspf', 'gmsppf'};

filterTypes  = {'ukf'};

date.day  = 17;
date.mon  = 11; 
date.year = 2017;

initialOrbit = loadInitialOrbit();
m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);

tStart = '00:00:00.000';
tEnd = '00:50:00.000';

timeData = TimeExt(tStart, tEnd, 10, date, 1e12); % time data for integrated navigation system
timeDataSubSystem  = TimeExt(tStart, tEnd, 1, date, 1e12); % time data for inertial navigation system & satellite navigation system

iterationNumber             = 1;
drawIterations              = 1;
esitimatedParams            = 6;

mass                        = 16;
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 1e-8*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 1e-6*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-7*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 1e-6*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-5*eye(3);      % [-]
levelArm                    = zeros(3, 1);      % [g] where g is Earth gravity acc
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);

initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
initialQuaternion           = initialOrbit(7:10);   % [-]

accelerationSigma           = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-7*ones(3, 1);      % [rad/sec]

insTrajInitErrorKm          = 5*ones(3, 1);         % [km]
insVelInitErrorKmSec        = 5e-3*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 1e-9*ones(4, 1);      % [-]

accScaleFactorInitError     = 1e-5*ones(3, 1);      % [-]
gyroScaleFactorInitError    = 1e-5*ones(3, 1);      % [-]
gyroBiasSigmaInitError      = 1e-3*gyroBiasSigma;   % [rad / sec]
accBiasSigmaInitError       = 1e-3*accBiasSigma;    % [km / sec^2]
reconciliationTime          = 1e12;                 % [sec]

ansSigmaTrajectory          = (1.1)^2;              % [km^2]  (3.45)^2   (1.1)^2
ansMeanTrajectory           = 0;                    % [km]
ansSigmaVelocity            = (1e-3)^2;             % [(km / sec)^2] (1.42e-3)^2  (1e-3)^2
ansMeanVelocity             = 0;                    % [(km / sec)]

errors = zeros(length(filterTypes), esitimatedParams, timeDataSubSystem.SimulationNumber);

alpha   = (1)^2;
beta    = (1)^2;
gamma   = (1)^2;

tic;
for l = 1:length(filterTypes)
    estimatorType = filterTypes(l);
    
    initCov = [diag(insTrajInitErrorKm).^2 zeros(3, 19); ... distance error [km]^2
        zeros(3, 3) diag(insVelInitErrorKmSec).^2 zeros(3, 16); ... velocity error [km/sec]^26
        zeros(4, 6) diag(insQuaternionInitError).^2 zeros(4, 12); ... quaternion error [-] 2.--- 5e-3
        zeros(3, 10) diag(accBiasSigmaInitError).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2 --- 5e-2
        zeros(3, 13) diag(gyroBiasSigmaInitError).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
        zeros(3, 16) diag(accScaleFactorInitError).^2 zeros(3, 3); ... acceler scale factor [-]
        zeros(3, 19) diag(gyroScaleFactorInitError).^2 ... gyro scale factor [-]
        ];
    
    iinsProcCov = getIntegratedInsProcCov(estimatorType, accBiasSigma, gyroBiasSigma);
    %========================================================================================================
    
    iinsInitArgs.initialParams(1:3)           = accBiasMu;
    iinsInitArgs.initialParams(4:6)           = accBiasSigma;
    iinsInitArgs.initialParams(7:9)           = gyroBiasMu;
    iinsInitArgs.initialParams(10:12)         = gyroBiasSigma;
    iinsInitArgs.initialParams(13:15)         = initialAcceleration;
    iinsInitArgs.initialParams(16:18)         = initialAngularVelocity;
    iinsInitArgs.initialParams(19:28)         = loadInitialOrbit();
    iinsInitArgs.initialParams(29)            = timeData.SampleTime;
    iinsInitArgs.initialParams(30)            = timeData.StartSecond;
    iinsInitArgs.processNoiseMean             = zeros(22, 1);
    iinsInitArgs.processNoiseMean(11:16)      = [accBiasMu; gyroBiasMu];
    iinsInitArgs.processNoiseCovariance       = beta*iinsProcCov;
    iinsInitArgs.observationNoiseMean         = [zeros(3, 1); zeros(3, 1)];
    iinsInitArgs.observationNoiseCovariance   = gamma*[[ansSigmaTrajectory*eye(3) zeros(3)]; [zeros(3) ansSigmaVelocity*eye(3)]];
    
    iterations = zeros(iterationNumber, esitimatedParams, timeDataSubSystem.SimulationNumber);
    x_iterations = zeros(iterationNumber, 10, timeDataSubSystem.SimulationNumber);
    
    fprintf('estimator: %s\n', estimatorType{1});
    
    gravModel = m_fitSolarSystemGravityModel(timeDataSubSystem.SampleTime, timeDataSubSystem.SimulationNumber);
    % parfor j = 1:iterationNumber
    for j = 1:iterationNumber
        insInitialState = initialOrbit + [insTrajInitErrorKm.*randn(3, 1); insVelInitErrorKmSec.*randn(3, 1); [1; 0; 0; 0] + insQuaternionInitError.*randn(4, 1)];
        insInitialState(7:10) = quaternionNormalize(insInitialState(7:10));
        
        accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataSubSystem, initialAcceleration, accelerationSigma);
        angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataSubSystem, initialAngularVelocity, angularVelocitySigma);
        
        simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataSubSystem, mass);
        xTrueState = simulator.simulate(iterationNumber == 1);
        insInitArgs = struct('accBiasMu', accBiasMu, 'accBiasSigma', accBiasSigma, 'gyroBiasMu', gyroBiasMu, 'gyroBiasSigma', gyroBiasSigma, 'visualize', 0,...
            'timeData', timeDataSubSystem, 'accNoiseVar', accNoiseVar, 'gyroNoiseVar', gyroNoiseVar, 'accScale', accScale, 'gyroScale', gyroScale,...
            'levelArm', levelArm, 'angularAccelerBodyFrame', angularAccelerBodyFrame, 'gyroGSensitiveBias', gyroGSensitiveBias,...
            'gravityModel', gravModel,'mass', mass, 'accelerationInBodyFrame', accelerationInBodyFrame, 'angularVelocityInBodyFrame', angularVelocityInBodyFrame);
        
        trueState = xTrueState;
        insSnsInitState = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]; [0; 0; 0]; 5e-6*ones(3, 1); 5e-6*ones(3, 1)];
        
        % gaussian | markov1 | wiener | exp | realImitator
        fns = NavImitator(trueState.State, ansSigmaTrajectory, ansMeanTrajectory, ansSigmaVelocity, ansMeanVelocity, 'gaussian');
        ins = initInertialNavigationSystem('init', insInitArgs);
        
        iins = IntegratedIns(ins, fns, timeData, iinsInitArgs, timeDataSubSystem, reconciliationTime);
        iinsState = iins.evaluate(insSnsInitState, alpha*initCov, insInitialState, estimatorType, iterationNumber == 1, iterationNumber == 1);
        
        quatErr = trueState.Rotation - iinsState.Rotation;
        errTraj = vectNormError(trueState.Trajectory, iinsState.Trajectory, 1e3);
        errVel  = vectNormError(trueState.Velocity, iinsState.Velocity, 1e3);
        
        iterations(j, :, :) = [errTraj; errVel; quatErr];
        x_iterations(j, :, :) = xTrueState.FullState - iinsState.FullState;
        
        fprintf('iteration of %d: completed\n', j);
    end
    
    if iterationNumber > 1
        for i = 1:timeDataSubSystem.SimulationNumber
            errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
        end
    end
    
    if drawIterations
        figure();
        for j = 1:min(iterationNumber, 200) % draw only 200 iterations due to laptop graphic card restrictions
            x_est = squeeze(x_iterations(j, :, :));
            
            subplot(1, 2, 1);
            plot2(timeDataSubSystem.TimeInMinutes, x_est(1:3, :)*1e3, sprintf('Ensemble of trajectory estimation errors (%s)', estimatorType{1}), ...
                {'Trajectory error'}, 'r error, meter', 'time, minutes');
            hold on;
            
            subplot(1, 2, 2);
            plot2(timeDataSubSystem.TimeInMinutes, x_est(4:6, :)*1e3, sprintf('Ensemble of velocity estimation errors (%s)', estimatorType{1}), ...
                {'Velocity error'}, 'v error, meter / sec', 'time, minutes');
            hold on;
        end
        
        figure();
        for j = 1:min(iterationNumber, 200)
            x_est = squeeze(x_iterations(j, :, :));
            plot2(timeDataSubSystem.TimeInMinutes, x_est(7:10, :), sprintf('Ensemble of quaternion estimation errors (%s)', estimatorType{1}), ...
                {'Quaternion error'}, 'q error', 'time, minutes');
            hold on;
        end
        
        figure();
        subplot(2, 1, 1);
        e1 = x_est(1:3, :)*1e3;
        plot2(timeDataSubSystem.TimeInMinutes, e1, sprintf('Transient response of estimation of trajectory error (%s)', estimatorType{1}),...
            {'r_x', 'r_y', 'r_z'}, 'r error, m', 'time, minutes');
        subplot(2, 1, 2);
        e2 = x_est(4:6, :)*1e3;
        plot2(timeDataSubSystem.TimeInMinutes, e2, sprintf('Transient response of estimation of velocity error (%s)', estimatorType{1}), ...
            {'v_x', 'v_y', 'v_z'}, 'v error, m / s', 'time, minutes');
    end
end

toc;

if iterationNumber > 1
    figure();
    subplot(2, 1, 1);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 1, :)), 'RMS trajectory', filterTypes, 'RMS_r, m', 'time, minutes');
    subplot(2, 1, 2);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 2, :)), 'RMS velocity', filterTypes, 'RMS_v, m / s', 'time, minutes');
    
    figure();
    subplot(4, 1, 1);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 3, :)), 'Quaternion error', filterTypes, 'q_0 error', 'time, minutes');
    subplot(4, 1, 2);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 4, :)), 'Quaternion error', filterTypes, 'q_1 error', 'time, minutes');
    subplot(4, 1, 3);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 5, :)), 'Quaternion error', filterTypes, 'q_2 error', 'time, minutes');
    subplot(4, 1, 4);
    plot2(timeDataSubSystem.TimeInMinutes, squeeze(errors(:, 6, :)), 'Quaternion error', filterTypes, 'q_3 error', 'time, minutes');
    
    for i = 1:length(filterTypes)
        estimatorType = filterTypes(i);
        fprintf('estimator: %s\n', estimatorType{1});
        fprintf('RMS trajectory: %d\n', errors(i, 1, end));
        fprintf('RMS velocity: %d\n', errors(i, 2, end));
        fprintf('RMS q_0: %d\n', errors(i, 3, end));
        fprintf('RMS q_1: %d\n', errors(i, 4, end));
        fprintf('RMS q_2: %d\n', errors(i, 5, end));
        fprintf('RMS q_3: %d\n', errors(i, 6, end));
    end
end
