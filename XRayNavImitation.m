close all force; clc;

addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('TOAMeasurementUnit')); % include folder with Pulsar & Quasar navigation
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));
addpath(genpath('Ephemeris'));
addpath(genpath('XNAV'));

secondInOneMinute       = 60;
sampleTime              = 0.0001; % seconds
simulationNumber        = 1*0.01/sampleTime; % x minutes * y seconds / sampleTime
simulationTime          = sampleTime*simulationNumber;
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / secondInOneMinute;
iterationNumber         = 1;

filterTypeArray         = {'ckf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf'};

T_till_current_epoch = 0.1465;
xRaySourceCount      = 7;
backgroundPhotnRate  = 5.9e4;
timeBucket           = 1e5; % 1e5
detectorArea         = 1;

earthEphemeris = loadEphemeris('earth', simulationNumber, secondInOneMinute/sampleTime);
sunEphemeris   = loadEphemeris('sun', simulationNumber, secondInOneMinute/sampleTime);
xRaySources    = loadXRaySources(xRaySourceCount);

initial = loadInitialOrbit();
initial = initial(1:6);

% simulate real trajectory of spaceship
simulator             = TrajectoryPhaseSpaceSatelliteSimulatorFree(initial, T_till_current_epoch);
spaceshipStateTrue    = simulator.Simulate(time);
initialSpaceshipState = spaceshipStateTrue(1:6, 1);

initialCov = [(1)^2*eye(3), zeros(3, 3); ...
    zeros(3, 3), (5e-3)^2*eye(3);
];
initialState = initial + chol(initialCov)*randn(6, 1);

estimatorType = filterTypeArray(1);

xRayDetector  = XRayDetector(xRaySources, detectorArea, timeBucket, backgroundPhotnRate, initial, earthEphemeris, sunEphemeris, T_till_current_epoch);
xRayNavSystem = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, xRayDetector);
tic;
stateVector = xRayNavSystem.Simulate(initialState, initialCov, time, estimatorType, T_till_current_epoch, sampleTime, 1);
fprintf('xRayNavSystem.Simulate: ')
toc;

figure();
    loglog(timeMinutes', abs(1e3*(stateVector(1:3, :) - spaceshipStateTrue(1:3, :))));
    ylabel('distance error, m');
    xlabel('time, minutes')
    title('distance error');
    legend('x axis', 'y axis', 'z axis');
    grid on;

figure();
    loglog(timeMinutes', abs(1e3*(stateVector(4:6, :) - spaceshipStateTrue(4:6, :))));
    ylabel('velocity error, m / sec')
    xlabel('time, minutes')
    title('velocity error');
    legend('x axis', 'y axis', 'z axis');
    grid on;