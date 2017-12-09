close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('../../'));

date.day  = 17;
date.mon  = 11;
date.year = 2017;

m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);

timeData = TimeExt('00:00:00.000', '01:00:00.000', 1, date, 1e9);
iterationNumber             = 1;
mass                        = 200; % [kg]

% large angular errors
%%{
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 5e-7*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 1e-4*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-3*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 1e-2*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-3*eye(3);      % [-]
levelArm                    = zeros(3, 1);
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);
initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
accelerationSigma           = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-2*ones(3, 1);      % [rad/sec]
insTrajInitErrorKm          = 3e-2*ones(3, 1);      % [km]
insVelInitErrorKmSec        = 5e-5*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 1e-3*ones(4, 1);      % [-]
accelerationInBodyFrame     = AccelerationInBodyFrame(timeData, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeData, initialAngularVelocity, angularVelocitySigma);
%}

% small angular errors
%{
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 5e-8*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 1e-5*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-6*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 1e-4*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-5*eye(3);      % [-]
levelArm                    = zeros(3, 1);
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);
initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
accelerationSigma           = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-4*ones(3, 16);      % [rad/sec]
insTrajInitErrorKm          = 3e-2*ones(3, 1);      % [km]
insVelInitErrorKmSec        = 5e-5*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 1e-3*ones(4, 1);      % [-]
accelerationInBodyFrame     = AccelerationInBodyFrame(timeData, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeData, initialAngularVelocity, angularVelocitySigma);
%}

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.visualize                    = 1;
insInitArgs.timeData                     = timeData;
insInitArgs.accNoiseVar                  = accNoiseVar;
insInitArgs.gyroNoiseVar                 = gyroNoiseVar;
insInitArgs.accScale                     = accScale;
insInitArgs.gyroScale                    = gyroScale;
insInitArgs.levelArm                     = levelArm;
insInitArgs.angularAccelerBodyFrame      = angularAccelerBodyFrame;
insInitArgs.gyroGSensitiveBias           = gyroGSensitiveBias;
insInitArgs.gravityModel                 = m_fitSolarSystemGravityModel(timeData.SampleTime, timeData.SimulationNumber);
insInitArgs.mass                         = mass;

iterations      = zeros(iterationNumber, 10, timeData.SimulationNumber);
iterationsErr   = zeros(iterationNumber, 6, timeData.SimulationNumber);
trueStateMatrix = zeros(iterationNumber, 10, timeData.SimulationNumber);

figure();
for i = 1:iterationNumber
    initialOrbit = loadInitialOrbit();
    initialOrbit(1:3) = initialOrbit(1:3) + 1e-1*randn(3, 1);
    initialOrbit(4:6) = initialOrbit(4:6) + 1e-3*randn(3, 1);
    initialOrbit(7:10) = initialOrbit(7:10) + 1e-5*randn(4, 1);
    initialOrbit(7:10) = quaternionNormalize(initialOrbit(7:10));
    
    simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeData, mass);
    trueState = simulator.simulate();
    
    insInitialState = initialOrbit + [insTrajInitErrorKm.*randn(3, 1); insVelInitErrorKmSec.*randn(3, 1); [1; 0; 0; 0] + insQuaternionInitError.*randn(4, 1)];
    insInitialState(7:10) = quaternionNormalize(insInitialState(7:10));
    
    ins = initInertialNavigationSystem('init', insInitArgs);
    
    estimations = ins.simulate(insInitialState);
    return
    iterations(i, :, :) = estimations;
    trueStateMatrix(i, :, :) = trueState.FullState;
    
    errTraj = vectNormError(trueState.Trajectory, estimations(1:3, :), 1e3);
    errVel  = vectNormError(trueState.Velocity, estimations(4:6, :), 1e3);
    errQuat = trueState.Rotation - estimations(7:10, :);
    iterationsErr(i, :, :) = [errTraj; errVel; errQuat];
    
    subplot(3, 1, 1);
    e1 = (trueState.Trajectory - estimations(1:3, :))*1e3;
    plot2(timeData.TimeInHour, e1, 'Ensemble of trajectory estimation errors', {'Trajectory error'}, 'Trajectory error, meter', 'time, hours');
    hold on;
    
    subplot(3, 1, 2);
    e2 = (trueState.Velocity - estimations(4:6, :))*1e3;
    plot2(timeData.TimeInHour, e2, 'Ensemble of velocity estimation errors', {'Velocity error'}, 'Velocity error, meter / sec', 'time, hours');
    hold on;
    
    subplot(3, 1, 3);
    plot2(timeData.TimeInHour, errQuat, 'Ensemble of angle rotation error', {'q_0', 'q_1', 'q_2', 'q_3'}, 'Quaternion error', 'time, hours');
    hold on;
    
    disp(i);
end

if iterationNumber > 1
    errors = zeros(6, timeData.SimulationNumber);
    
    for i = 1:timeData.SimulationNumber
        errors(:, i) = ( sum( ( squeeze(iterationsErr(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    
    rErrors = squeezeToRow(errors(1, :));
    vErrors = squeezeToRow(errors(2, :));
    
    figure();
    subplot(2, 1, 1);
    plot2(timeData.TimeInHour, rErrors, 'RMS trajectory', {'INS'}, 'RMS trajectory error, meter');
    subplot(2, 1, 2);
    plot2(timeData.TimeInHour, vErrors, 'RMS velocity', {'INS'}, 'RMS velocity error, meter / sec');
    
    figure();
    subplot(4, 1, 1);
    plot2(timeData.TimeInHour, squeeze(errors(3, :)), 'RMS quaternion error', {'INS'}, 'q_0');
    subplot(4, 1, 2);
    plot2(timeData.TimeInHour, squeeze(errors(4, :)), 'RMS quaternion error', {'INS'}, 'q_1');
    subplot(4, 1, 3);
    plot2(timeData.TimeInHour, squeeze(errors(5, :)), 'RMS quaternion error', {'INS'}, 'q_2');
    subplot(4, 1, 4);
    plot2(timeData.TimeInHour, squeeze(errors(6, :)), 'RMS quaternion error', {'INS'}, 'q_3');
end

pdfEvolutionPlot2(iterations(:, 1, :), trueStateMatrix(:, 1, :), {'r_X'}, {', m'}, timeData);
pdfEvolutionPlot2(iterations(:, 2, :), trueStateMatrix(:, 2, :), {'r_Y'}, {', m'}, timeData);
pdfEvolutionPlot2(iterations(:, 3, :), trueStateMatrix(:, 3, :), {'r_Z'}, {', m'}, timeData);

pdfEvolutionPlot2(iterations(:, 4, :), trueStateMatrix(:, 4, :), {'v_X'}, {', m / s'}, timeData);
pdfEvolutionPlot2(iterations(:, 5, :), trueStateMatrix(:, 5, :), {'v_Y'}, {', m / s'}, timeData);
pdfEvolutionPlot2(iterations(:, 6, :), trueStateMatrix(:, 6, :), {'v_Z'}, {', m / s'}, timeData);

pdfEvolutionPlot2(iterations(:, 7, :), trueStateMatrix(:, 7, :), {'q_0'}, {''}, timeData);
pdfEvolutionPlot2(iterations(:, 8, :), trueStateMatrix(:, 8, :), {'q_1'}, {''}, timeData);
pdfEvolutionPlot2(iterations(:, 9, :), trueStateMatrix(:, 9, :), {'q_2'}, {''}, timeData);
pdfEvolutionPlot2(iterations(:, 10, :), trueStateMatrix(:, 10, :), {'q_3'}, {''}, timeData);

pdfPlot2(iterations(:, 1:3, :), trueStateMatrix(:, 1:3, :), {'r_X', 'r_Y', 'r_Z'});
pdfPlot2(iterations(:, 4:6, :), trueStateMatrix(:, 4:6, :), {'v_X', 'v_Y', 'v_Z'});
pdfPlot2(iterations(:, 7:10, :), trueStateMatrix(:, 7:10, :), {'q_0', 'q_1', 'q_2', 'q_2'});
