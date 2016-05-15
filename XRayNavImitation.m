close all force; clc;

addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('TOAMeasurementUnit')); % include folder with Pulsar & Quasar navigation
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));
addpath(genpath('Ephemeris'));
addpath(genpath('XNAV'));

sampleTime              = 1; % seconds
simulationNumber        = 1e3;
simulationTime          = sampleTime*simulationNumber;
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;

filterTypeArray         = {'ckf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf'};

T_till_current_epoch = 0.1465;
xRaySourceCount      = 7;
backgroundPhotnRate  = 5.9;
timeBucket           = 1e-1; % 1e5
detectorArea         = 1;

earthEphemeris = loadEphemeris('earth', simulationNumber);
sunEphemeris = loadEphemeris('sun', simulationNumber);

xRaySources = loadXRaySources(xRaySourceCount);
initialSpaceshipState = loadInitialOrbit();
initialSpaceshipState = initialSpaceshipState(1:6);

% simulate real trajectory of spaceship
simulator = TrajectoryPhaseSpaceSatelliteSimulatorFree(initialSpaceshipState, T_till_current_epoch);
spaceshipStateTrue = simulator.Simulate(time);

xRayDetector = XRayDetector(xRaySources, detectorArea, timeBucket, backgroundPhotnRate, initialSpaceshipState, earthEphemeris, sunEphemeris, T_till_current_epoch);

initialCov = [(100)^2*eye(3), zeros(3, 3); ...
    zeros(3, 3), (5e-1)^2*eye(3);
];
initialState = spaceshipStateTrue(:, 1) + chol(initialCov)*randn(6, 1);

estimatorType = filterTypeArray(1);

xRayNavSystem = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, xRayDetector);
stateVector = xRayNavSystem.Simulate(initialState, initialCov, time, estimatorType, T_till_current_epoch, sampleTime);

figure();
    plot(timeMinutes', stateVector(1:3, :) - spaceshipStateTrue(1:3, :));
    title('distance error');
    legend('x axis', 'y axis', 'z axis');
    grid on;

figure();
    plot(timeMinutes', stateVector(4:6, :) - spaceshipStateTrue(4:6, :));
    title('velocity error');
    legend('x axis', 'y axis', 'z axis');
    grid on;