close all force; clc; clearvars;

addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('InertialMeasurementUnit')); % include folder with inertial measurement functions
addpath(genpath('InertialNavigationSystem')); % include folder with inertial navigation system functions
addpath(genpath('SNS')); % include folder with inertial navigation system functions
addpath(genpath('Utils'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;

initialOrbit = loadInitialOrbit();

timeData = TimeExt('00:00:00.000', '05:00:00.000', 1e-1, date, 1);
iterationNumber             = 200;
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 5e-6*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 5e-4*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-5*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 5e-4*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-5*eye(3);      % [-]
levelArm                    = zeros(3, 1);
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);
initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
initialQuaternion           = initialOrbit(7:10);   % [-]
accelerationSigma           = 2e-4*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-3*ones(3, 1);      % [rad/sec]
insTrajInitErrorKm          = 3e-1*ones(3, 1);      % [km]
insVelInitErrorKmSec        = 5e-4*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 5e-5*ones(4, 1);      % [-] error converted to angle approx equal 3-5 grad.
accScaleFactorInitError     = 1e-5*ones(3, 1);      % [-]
gyroScaleFactorInitError    = 1e-5*ones(3, 1);      % [-]
gyroBiasSigmaInitError      = 1e-2*gyroBiasSigma;   % [rad / sec]
accBiasSigmaInitError       = 1e-2*accBiasSigma;    % [km / sec^2]

accelerationInBodyFrame     = AccelerationInBodyFrame(timeData, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeData, initialAngularVelocity, angularVelocitySigma);


simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeData);
starshipTrueState = simulator.simulate();

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.visualize                    = 0;
insInitArgs.timeData                     = timeData;
insInitArgs.accNoiseVar                  = accNoiseVar;
insInitArgs.gyroNoiseVar                 = gyroNoiseVar;
insInitArgs.accScale                     = accScale;
insInitArgs.gyroScale                    = gyroScale;
insInitArgs.levelArm                     = levelArm;
insInitArgs.angularAccelerBodyFrame      = angularAccelerBodyFrame;
insInitArgs.gyroGSensitiveBias           = gyroGSensitiveBias;

iterations = zeros(iterationNumber, 10, timeData.SimulationNumber);

parfor i = 1:iterationNumber
    insInitialState = initialOrbit + [insTrajInitErrorKm.*randn(3, 1); insVelInitErrorKmSec.*randn(3, 1); [1; 0; 0; 0] + insQuaternionInitError.*randn(4, 1)];
    insInitialState(7:10) = quaternionNormalize(insInitialState(7:10));
    
    ins = initInertialNavigationSystem('init', insInitArgs);
    iterations(i, :, :) = ins.simulate(insInitialState);
    disp(i);
end

legend = {'r_X', 'r_Y', 'r_Z', 'v_X', 'v_Y', 'v_Z', 'q_0', 'q_X', 'q_Y', 'q_Z'};
for i = 1:10
    figure();
    sample = squeeze(iterations(:, i, end)) - rvecrep(starshipTrueState.FullState(i, end), iterationNumber);
%     hist(sample, 'Normalization', 'pdf' );
    h1 = histfit(sample, 20, 'kernel');
    h1(1).FaceColor = [.8 .8 1];
    h1(2).Color = [.2 .2 .2];
    
    hold on
    h2 = histfit(sample, 20, 'normal');
    
    h2(2).Color = [.7 .7 .7];
    delete(h2(1))
    
    title(strcat(legend{i}, ' pdf'));
end


