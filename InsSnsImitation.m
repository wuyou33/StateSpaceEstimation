close all force; clc; clearvars;
addpath(genpath('./'));

% ghqf not working due to dimension, too number of points :(
%{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf', 'srcdkf', 'pf', 'sppf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'gspf', 'gmsppf'};
estimatorType  = {'gmsppf'};

date.day  = 17;
date.mon  = 11;
date.year = 2015;

initialOrbit = loadInitialOrbit();

tStart = '00:00:00.000';
% tEnd = '00:00:20.000';
tEnd = '00:00:00.200';

timeData = TimeExt(tStart, tEnd, 1e-3, date, 1); % time data for integrated navigation system
timeDataSubSystem  = TimeExt(tStart, tEnd, 1e-4, date, 1); % time data for inertial navigation system & satellite navigation system

iterationNumber             = 1;
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
insQuaternionInitError      = 5e-5*ones(4, 1);      % [-] error converted to angle approx equal 3-5 grad.
accScaleFactorInitError     = 1e-5*ones(3, 1);      % [-]
gyroScaleFactorInitError    = 1e-5*ones(3, 1);      % [-]
gyroBiasSigmaInitError      = 1e-2*gyroBiasSigma;   % [rad / sec]
accBiasSigmaInitError       = 1e-2*accBiasSigma;    % [km / sec^2]

accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataSubSystem, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataSubSystem, initialAngularVelocity, angularVelocitySigma);

simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataSubSystem);
starshipTrueState = simulator.simulate();

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.visualize                    = 1;
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

%========================================================================================================
initCov = [diag(insTrajInitErrorKm).^2 zeros(3, 19); ... distance error [km]^2
    zeros(3, 3) diag(insVelInitErrorKmSec).^2 zeros(3, 16); ... velocity error [km/sec]^2
    zeros(4, 6) diag(insQuaternionInitError).^2 zeros(4, 12); ... quaternion error [-] 2.--- 5e-3
    zeros(3, 10) diag(accBiasSigmaInitError).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2 --- 5e-2
    zeros(3, 13) diag(gyroBiasSigmaInitError).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
    zeros(3, 16) diag(1e-3*accScaleFactorInitError).^2 zeros(3, 3); ... acceler scale factor [-]
    zeros(3, 19) diag(1e-3*gyroScaleFactorInitError).^2 ... gyro scale factor [-]
    ];

%========================================================================================================

% % nice setup for sigma point filter and gaussian noise.
% insSnsProcCov = [(1e-4*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
%                  zeros(3, 3) (5e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
%                  zeros(4, 6) (1e-7*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
%                  zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
%                  zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
%                  zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
%                  zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
%                 ];

% nice setup for particle filter and gaussian noise.
% insSnsProcCov = [(7.25e-4*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
%                 zeros(3, 3) (7.85e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
%                 zeros(4, 6) (5.75e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
%                 zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
%                 zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
%                 zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
%                 zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
%             ];

% nice setup for high(fifth)-degree ckf
% insSnsProcCov = [(1.25e-2*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
%                 zeros(3, 3) (7.85e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
%                 zeros(4, 6) (5.75e-5*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
%                 zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
%                 zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
%                 zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
%                 zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
%             ];

% nice setup for cubature-quadrature kalman filter and gaussian noise.
% insSnsProcCov = [(5.75e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
%                  zeros(3, 3) (5e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
%                  zeros(4, 6) (5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
%                  zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
%                  zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
%                  zeros(3, 16) (1e-8*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
%                  zeros(3, 19) (1e-8*eye(3)).^2 ... gyro scale factor [-]
%                 ];

% nice setup for gaussian sum particle filter and gaussian noise.
insSnsProcCov = [(1e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
                 zeros(3, 3) (1e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
                 zeros(4, 6) (2.5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
                 zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
                 zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
                 zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
                 zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
                ];

% nice setup for gaussian sum sigma point particle filter and gaussian noise.
% insSnsProcCov = [(1e-4*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
%                  zeros(3, 3) (2.5e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
%                  zeros(4, 6) (2.5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
%                  zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
%                  zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
%                  zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
%                  zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
%                 ];
            
rCovKoeff = 1.0; % 1.5 for pf, otherwise 1
vCovKoeff = 1.0; % 1.5 for pf, otherwise 1
insSnsInitArgs.initialParams(1:3)           = accBiasMu;
insSnsInitArgs.initialParams(4:6)           = accBiasSigma;
insSnsInitArgs.initialParams(7:9)           = gyroBiasMu;
insSnsInitArgs.initialParams(10:12)         = gyroBiasSigma;
insSnsInitArgs.initialParams(13:15)         = initialAcceleration;
insSnsInitArgs.initialParams(16:18)         = initialAngularVelocity;
insSnsInitArgs.initialParams(19:22)         = initialQuaternion;
insSnsInitArgs.initialParams(23)            = timeData.SampleTime;
insSnsInitArgs.initialParams(24)            = timeData.StartSecond;
insSnsInitArgs.processNoiseMean             = zeros(22, 1);
insSnsInitArgs.processNoiseMean(11:16)      = [accBiasMu; gyroBiasMu];
insSnsInitArgs.processNoiseCovariance       = insSnsProcCov;
insSnsInitArgs.observationNoiseMean         = [zeros(3, 1); zeros(3, 1)]; % [1e-2*ones(3, 1); zeros(3, 1)];
insSnsInitArgs.observationNoiseCovariance   = [rCovKoeff*[(1e-2*eye(3)).^2 zeros(3)]; vCovKoeff*[zeros(3) (1e-5*eye(3)).^2]];

esitimatedParams = 5;
iterations = zeros(iterationNumber, esitimatedParams, timeDataSubSystem.SimulationNumber);

fprintf('estimator: %s\n', estimatorType{1});
realTime = timeDataSubSystem.Time;
%     parfor j = 1:iterationNumber
for j = 1:iterationNumber
    insSnsInitState = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]; [0; 0; 0]; 5e-5*ones(3, 1); 5e-5*ones(3, 1)];
    
    sns = SnsImitator(starshipTrueState.State, 'gaussian'); % gaussian markov1 wiener exp
    ins = initInertialNavigationSystem('init', insInitArgs);
    
    insSns = IntegratedInsSns(ins, sns, timeData, insSnsInitArgs, timeDataSubSystem);
    tic;
    insSnsState = insSns.evaluate(insSnsInitState, initCov, insInitialState, estimatorType, 0, 1);
    toc;
    angErr = angleErrorsFromQuaternion(insSnsState.Rotation, starshipTrueState.Rotation);
    errTraj = [vectNormError(starshipTrueState.Trajectory, sns.Trajectory, 1e3); ...
        vectNormError(starshipTrueState.Trajectory, insSnsState.Trajectory, 1e3)];
    errVel  = [vectNormError(starshipTrueState.Velocity, sns.Velocity, 1e3); ...
        vectNormError(starshipTrueState.Velocity, insSnsState.Velocity, 1e3)];
    
    iterations(j, :, :) = [errTraj(2, :); errVel(2, :); angErr];
        
    if iterationNumber == 1 && length(estimatorType) == 1
        figure();
        subplot(2, 1, 1);
        e1 = (starshipTrueState.Trajectory - insSnsState.Trajectory)*1e3;
        plot2(realTime, e1, 'trajectory INS & SNS error', {'x', 'y', 'z'}, 'trajectory error, m');
        subplot(2, 1, 2);
        e2 = (starshipTrueState.Velocity - insSnsState.Velocity)*1e3;
        plot2(realTime, e2, 'velocity INS & SNS error', {'x', 'y', 'z'}, 'velocity error, m/sec');
        
        figure();
        subplot(2, 1, 1);
        plot2(realTime, errTraj, 'trajectory errors in SNS and SNS & INS', {'SNS', 'INS & SNS'}, 'trajectory error, meter');
        subplot(2, 1, 2);
        plot2(realTime, errVel, 'velocity errors in SNS and SNS & INS', {'SNS', 'INS & SNS'}, 'velocity error, meter / sec');
        
        figure()
        plot2(realTime, angErr, 'angle rotation error', {'yaw', 'pitch', 'roll'}, 'angle error, rad');
    end
    
    fprintf('iteration of %d: completed\n', j);
end

if iterationNumber > 1
    errors = zeros(esitimatedParams, timeData.SimulationNumber);
    for i = 1:timeData.SimulationNumber
        errors(:, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    
    figure();
    subplot(3, 1, 1);
    plot2(realTime, errors(1, :), 'trajectory errors in SNS & INS', {'INS & SNS'}, 'trajectory error, meter');
    subplot(3, 1, 2);
    plot2(realTime, errors(2, :), 'velocity errors in SNS & INS', {'INS & SNS'}, 'velocity error, meter / sec');
    subplot(3, 1, 3);
    plot2(realTime, errors(3:5, :), 'angle error and SNS & INS', {'yaw', 'pitch', 'roll'}, 'velocity error, meter / sec');
end