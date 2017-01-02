close all force; clc; clearvars;
addpath(genpath('./'));

% ghqf not working due to dimension, too number of points :(
% need to check point and weights generation in 'sghqf'
%{'ekf', 'ukf', 'cdkf', 'ckf', 'sckf', 'srukf', 'srcdkf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'pf', 'sppf', 'gspf', 'gmsppf'};

% non-linear function approximation filters
% filterTypes  = {'ekf'};

% square root sigma points family
% filterTypes  = {'sckf', 'srukf', 'srcdkf'};

% sigma points family
% filterTypes  = {'ckf', 'ukf', 'cdkf'};

% high degree sigma points family
% filterTypes  = {'fdckf', 'cqkf', 'sghqf'};

% particle filters family
% filterTypes  = {'pf', 'sppf', 'gspf', 'gmsppf'};

filterTypes  = {'srukf'};

date.day  = 17;
date.mon  = 11;
date.year = 2015;

initialOrbit = loadInitialOrbit();

tStart = '00:00:00.000';
tEnd = '00:05:00.000';

timeData = TimeExt(tStart, tEnd, 1e-1, date, 1e9); % time data for integrated navigation system 
timeDataSubSystem  = TimeExt(tStart, tEnd, 1e-1, date, 1e9); % time data for inertial navigation system & satellite navigation system

logLastError                = 1;
iterationNumber             = 1;
isOneLaunch                 = iterationNumber == 1;
esitimatedParams            = 6;
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 5e-8*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 1e-6*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-6*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 1e-5*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-5*eye(3);      % [-]
levelArm                    = zeros(3, 1);
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);
initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
initialQuaternion           = initialOrbit(7:10);   % [-]
accelerationSigma           = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-4*ones(3, 1);      % [rad/sec]
insTrajInitErrorKm          = 3e-2*ones(3, 1);      % [km]
insVelInitErrorKmSec        = 5e-5*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 1e-3*ones(4, 1);      % [-]
accScaleFactorInitError     = 1e-5*ones(3, 1);      % [-]
gyroScaleFactorInitError    = 1e-5*ones(3, 1);      % [-]
gyroBiasSigmaInitError      = 1e-2*gyroBiasSigma;   % [rad / sec]
accBiasSigmaInitError       = 1e-2*accBiasSigma;    % [km / sec^2]
snsSigmaTrajectory          = 3*(1e-2)^2;           % [km^2]
snsMeanTrajectory           = 0;                    % [km]
snsSigmaVelocity            = 3*(1e-5)^2;           % [(km / sec)^2]
snsMeanVelocity             = 0;                    % [(km / sec)]
reconciliationTime          = 100;                  % [sec]
accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataSubSystem, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataSubSystem, initialAngularVelocity, angularVelocitySigma);

simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataSubSystem);
starshipTrueState = simulator.simulate(isOneLaunch);

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.visualize                    = 0;
insInitArgs.timeData                     = timeDataSubSystem;
insInitArgs.accNoiseVar                  = accNoiseVar;
insInitArgs.gyroNoiseVar                 = gyroNoiseVar;
insInitArgs.accScale                     = accScale;
insInitArgs.gyroScale                    = gyroScale;
insInitArgs.levelArm                     = levelArm;
insInitArgs.angularAccelerBodyFrame      = angularAccelerBodyFrame;
insInitArgs.gyroGSensitiveBias           = gyroGSensitiveBias;

insInitialState = initialOrbit + [insTrajInitErrorKm.*randn(3, 1); insVelInitErrorKmSec.*randn(3, 1); [1; 0; 0; 0] + insQuaternionInitError.*randn(4, 1)];
insInitialState(7:10) = quaternionNormalize(insInitialState(7:10));

errors = zeros(length(filterTypes), esitimatedParams, timeDataSubSystem.SimulationNumber);
initialConditionStabilityCoeff = 1;

tic;

for l = 1:length(filterTypes)
    estimatorType = filterTypes(l);
    
    %========================================================================================================
    initCov = [diag(insTrajInitErrorKm).^2 zeros(3, 19); ... distance error [km]^2
        zeros(3, 3) diag(insVelInitErrorKmSec).^2 zeros(3, 16); ... velocity error [km/sec]^2
        zeros(4, 6) diag(insQuaternionInitError).^2 zeros(4, 12); ... quaternion error [-] 2.--- 5e-3
        zeros(3, 10) diag(accBiasSigmaInitError).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2 --- 5e-2
        zeros(3, 13) diag(gyroBiasSigmaInitError).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
        zeros(3, 16) diag(1e-3*accScaleFactorInitError).^2 zeros(3, 3); ... acceler scale factor [-]
        zeros(3, 19) diag(1e-3*gyroScaleFactorInitError).^2 ... gyro scale factor [-]
        ];
    
    if stringmatch(estimatorType, {'ukf', 'cdkf', 'srukf', 'srcdkf', 'ckf', 'sckf', 'ekf', 'kf'})
        insSnsProcCov = [(1.5e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (5e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (1e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) 1e2*diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) 1e0*diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-9*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-9*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'pf', 'gspf', 'gmsppf', 'sppf'})
        insSnsProcCov = [(1e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (1e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (2.5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'sghqf'})
        insSnsProcCov = [(1.25e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (7.85e-5*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (5.75e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-9*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-9*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'cqkf', 'fdckf'})
        insSnsProcCov = [(5.75e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (5e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-8*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-8*eye(3)).^2 ... gyro scale factor [-]
            ];
    else
        error('please specify process noise covariance matrix for every filtration algorithm');
    end
    %========================================================================================================
    
    insSnsInitArgs.initialParams(1:3)           = accBiasMu;
    insSnsInitArgs.initialParams(4:6)           = accBiasSigma;
    insSnsInitArgs.initialParams(7:9)           = gyroBiasMu;
    insSnsInitArgs.initialParams(10:12)         = gyroBiasSigma;
    insSnsInitArgs.initialParams(13:15)         = initialAcceleration;
    insSnsInitArgs.initialParams(16:18)         = initialAngularVelocity;
    insSnsInitArgs.initialParams(19:28)         = insInitialState;
    insSnsInitArgs.initialParams(29)            = timeData.SampleTime;
    insSnsInitArgs.initialParams(30)            = timeData.StartSecond;
    insSnsInitArgs.processNoiseMean             = zeros(22, 1);
    insSnsInitArgs.processNoiseMean(11:16)      = [accBiasMu; gyroBiasMu];
    insSnsInitArgs.processNoiseCovariance       = insSnsProcCov;
    insSnsInitArgs.observationNoiseMean         = [zeros(3, 1); zeros(3, 1)]; % [1e-2*ones(3, 1); zeros(3, 1)];
    insSnsInitArgs.observationNoiseCovariance   = [[(1e-2*eye(3)).^2 zeros(3)]; [zeros(3) (1e-5*eye(3)).^2]];
    
    iterations = zeros(iterationNumber, esitimatedParams, timeDataSubSystem.SimulationNumber);
    
    fprintf('estimator: %s\n', estimatorType{1});
    realTime = timeDataSubSystem.Time;
    snsState = starshipTrueState.State;
    trueRotation   = starshipTrueState.Rotation;
    trueTrajectory = starshipTrueState.Trajectory;
    trueVelocity   = starshipTrueState.Velocity;
    
    for j = 1:iterationNumber
        insSnsInitState = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]; [0; 0; 0]; 5e-5*ones(3, 1); 5e-5*ones(3, 1)];        
        
        % gaussian | markov1 | wiener | exp | realImitator
        sns = SnsImitator(snsState, snsSigmaTrajectory, snsMeanTrajectory, snsSigmaVelocity, snsMeanVelocity, 'gaussian');
        ins = initInertialNavigationSystem('init', insInitArgs);
        
        insSns = IntegratedInsSns(ins, sns, timeData, insSnsInitArgs, timeDataSubSystem, reconciliationTime);
        insSnsState = insSns.evaluate(insSnsInitState, initialConditionStabilityCoeff*initCov, insInitialState, estimatorType, isOneLaunch, isOneLaunch);
        
        quatErr = trueRotation - insSnsState.Rotation;
        errTraj = vectNormError(trueTrajectory, insSnsState.Trajectory, 1e3);
        errVel  = vectNormError(trueVelocity, insSnsState.Velocity, 1e3);
        
        iterations(j, :, :) = [errTraj; errVel; quatErr];
        
        if isOneLaunch
            figure();
            subplot(2, 1, 1);
            e1 = (trueTrajectory - insSnsState.Trajectory)*1e3;
            plot2(realTime, e1, 'trajectory INS & SNS error', {'x', 'y', 'z'}, 'trajectory error, m');
            subplot(2, 1, 2);
            e2 = (trueVelocity - insSnsState.Velocity)*1e3;
            plot2(realTime, e2, 'velocity INS & SNS error', {'x', 'y', 'z'}, 'velocity error, m/sec');
            
            figure();
            subplot(2, 1, 1);
            plot2(realTime, errTraj, 'trajectory errors in SNS and SNS & INS', estimatorType, 'trajectory error, meter');
            subplot(2, 1, 2);
            plot2(realTime, errVel, 'velocity errors in SNS and SNS & INS', estimatorType, 'velocity error, meter / sec');
            
            figure()
            plot2(realTime, quatErr, 'angle rotation error', {'q0', 'q1', 'q2', 'q3'}, 'quaternion error');
        end
        
        fprintf('iteration of %d: completed\n', j);
    end
    
    if iterationNumber > 1
        for i = 1:timeDataSubSystem.SimulationNumber
            errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
        end
    end
    
    if logLastError
        fprintf('estimator: %s\n', estimatorType{1});
        fprintf('RMS trajectory: %d\n', errors(l, 1, end));
        fprintf('RMS velocity: %d\n', errors(l, 2, end));
    end
end

toc;

if iterationNumber > 1
    %     fileName = strcat('errors_ins_sns_pf_family.mat');
    %     save(fileName, 'errors');
    
    legends  = mergeCells({'sns'}, filterTypes);
    
    rErrors = [(sqrt( snsSigmaTrajectory ))*1e3*ones(1, timeDataSubSystem.SimulationNumber); squeezeToRow(errors(:, 1, :))];
    vErrors = [(sqrt( snsSigmaVelocity ))*1e3*ones(1, timeDataSubSystem.SimulationNumber); squeezeToRow(errors(:, 2, :))];
    %     rErrors = squeeze(errors(:, 1, :));
    %     vErrors = squeeze(errors(:, 2, :));
    
    figure();
    subplot(2, 1, 1);
    plot2(realTime, rErrors, 'trajectory errors in SNS & INS', legends, 'trajectory error, meter');
    subplot(2, 1, 2);
    plot2(realTime, vErrors, 'velocity errors in SNS & INS', legends, 'velocity error, meter / sec');
    
    figure();
    subplot(4, 1, 1);
    plot2(realTime, squeeze(errors(:, 3, :)), 'quaternion error in SNS & INS', filterTypes, 'q0 error');
    subplot(4, 1, 2);
    plot2(realTime, squeeze(errors(:, 4, :)), 'quaternion error in SNS & INS', filterTypes, 'q1 error');
    subplot(4, 1, 3);
    plot2(realTime, squeeze(errors(:, 5, :)), 'quaternion error in SNS & INS', filterTypes, 'q2 error');
    subplot(4, 1, 4);
    plot2(realTime, squeeze(errors(:, 6, :)), 'quaternion error in SNS & INS', filterTypes, 'q3 error');
end