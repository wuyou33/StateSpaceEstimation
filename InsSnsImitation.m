close all force; clc; clearvars;
addpath(genpath('./'));

% ghqf not working due to dimension, too number of points :(
%{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf', 'srcdkf', 'pf', 'sppf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'gspf', 'gmsppf'};

% square root sigma points family
% filterTypes  = {'sckf', 'srukf', 'srcdkf'};
% legends  = {'sns', 'sckf', 'srukf', 'srcdkf'};

% sigma points family
% filterTypes  = {'ckf', 'ukf', 'cdkf'};
% legends  = {'sns', 'ckf', 'ukf', 'cdkf'};

% high degree sigma points family
% filterTypes  = {'fdckf', 'cqkf', 'sghqf'};
% legends  = {'sns', 'fdckf', 'cqkf', 'sghqf'};

% particle filters family
filterTypes  = {'pf', 'sppf', 'gspf', 'gmsppf'};
legends  = {'sns', 'pf', 'sppf', 'gspf', 'gmsppf'};

date.day  = 17;
date.mon  = 11;
date.year = 2015;

initialOrbit = loadInitialOrbit();

tStart = '00:00:00.000';
% tEnd = '00:00:20.000';
tEnd = '00:00:0.500';

timeData = TimeExt(tStart, tEnd, 1e-3, date, 1); % time data for integrated navigation system
timeDataSubSystem  = TimeExt(tStart, tEnd, 1e-4, date, 1); % time data for inertial navigation system & satellite navigation system

iterationNumber             = 25;
esitimatedParams            = 5;
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
    
    if stringmatch(estimatorType, {'ukf', 'cdkf', 'srukf', 'srcdkf', 'ckf', 'sckf'})
        insSnsProcCov = [(2e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (2e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) 1e1*diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) 1e1*diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-8*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-8*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'pf', 'gspf'})
        insSnsProcCov = [(1e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (1e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (2.5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'gmsppf', 'sppf'})
        insSnsProcCov = [(1e-4*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (2.5e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (2.5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif stringmatch(estimatorType, {'fdckf'})
        insSnsProcCov = [(1.25e-2*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (7.85e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (5.75e-5*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
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
    elseif stringmatch(estimatorType, {'cqkf'})
        insSnsProcCov = [(5.75e-3*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (5e-6*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (5e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-8*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-8*eye(3)).^2 ... gyro scale factor [-]
            ];
    end
    %========================================================================================================
    
    % nice setup for particle filter and gaussian noise.
    % insSnsProcCov = [(7.25e-4*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
    %                 zeros(3, 3) (7.85e-7*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
    %                 zeros(4, 6) (5.75e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
    %                 zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
    %                 zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
    %                 zeros(3, 16) (1e-7*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
    %                 zeros(3, 19) (1e-7*eye(3)).^2 ... gyro scale factor [-]
    %             ];
    
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
    insSnsInitArgs.observationNoiseCovariance   = [[(1e-2*eye(3)).^2 zeros(3)]; [zeros(3) (1e-5*eye(3)).^2]];
    
    iterations = zeros(iterationNumber, esitimatedParams, timeDataSubSystem.SimulationNumber);
    
    fprintf('estimator: %s\n', estimatorType{1});
    realTime = timeDataSubSystem.Time;
    snsState = starshipTrueState.State;
    trueRotation   = starshipTrueState.Rotation;
    trueTrajectory = starshipTrueState.Trajectory;
    trueVelocity   = starshipTrueState.Velocity;
    
%     parfor j = 1:iterationNumber
    for j = 1:iterationNumber
        insSnsInitState = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]; [0; 0; 0]; 5e-5*ones(3, 1); 5e-5*ones(3, 1)];
        
        sns = SnsImitator(snsState, 'gaussian'); % gaussian markov1 wiener exp
        ins = initInertialNavigationSystem('init', insInitArgs);
        
        insSns = IntegratedInsSns(ins, sns, timeData, insSnsInitArgs, timeDataSubSystem);
        %         tic;
        insSnsState = insSns.evaluate(insSnsInitState, initCov, insInitialState, estimatorType, iterationNumber == 1, iterationNumber == 1);
        %         toc;
        angErr = angleErrorsFromQuaternion(insSnsState.Rotation, trueRotation);
        errTraj = vectNormError(trueTrajectory, insSnsState.Trajectory, 1e3);
        errVel  = vectNormError(trueVelocity, insSnsState.Velocity, 1e3);
        
        iterations(j, :, :) = [errTraj; errVel; angErr];
        
        if iterationNumber == 1
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
        for i = 1:timeDataSubSystem.SimulationNumber
            errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
        end
    end
end

if iterationNumber > 1
    fileName = strcat('errors_ins_sns_pf_family.mat');
    save(fileName, 'errors');
    
    rErrors = [(sqrt( 3*(1e-2)^2 ))*1e3*ones(1, timeDataSubSystem.SimulationNumber); squeeze(errors(:, 1, :))];
    vErrors = [(sqrt( 3*(1e-5)^2 ))*1e3*ones(1, timeDataSubSystem.SimulationNumber); squeeze(errors(:, 2, :))];
    figure();
    subplot(2, 1, 1);
    plot2(realTime, rErrors, 'trajectory errors in SNS & INS', legends, 'trajectory error, meter');
    subplot(2, 1, 2);
    plot2(realTime, vErrors, 'velocity errors in SNS & INS', legends, 'velocity error, meter / sec');
    
    figure();
    subplot(3, 1, 1);
    plot2(realTime, squeeze(errors(:, 3, :)), 'angle error and SNS & INS', filterTypes, 'yaw error, rad');
    subplot(3, 1, 2);
    plot2(realTime, squeeze(errors(:, 4, :)), 'angle error and SNS & INS', filterTypes, 'pitch error, rad');
    subplot(3, 1, 3);
    plot2(realTime, squeeze(errors(:, 5, :)), 'angle error and SNS & INS', filterTypes, 'roll error, rad');
end
