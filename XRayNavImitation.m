close all; clc; clearvars;

addpath(genpath('OrbitalMotion'));
addpath(genpath('Statistics'));
addpath(genpath('TOAMeasurementUnit'));
addpath(genpath('StateSpaceEstimation'));
addpath(genpath('Utils'));
addpath(genpath('Ephemeris'));
addpath(genpath('XNAV'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd = '00:00:00.500';
timeDataXRay  = TimeExt(timeStart, timeEnd, 1e-3, date, 1);
iterationNumber    = 1;
secondInOneMinute  = 60;

filterType = {'srukf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf', 'sppf', 'fdckf', 'fdckfAugmented', 'cqkf', 'gspf', 'gmsppf', 'ghqf', 'sghqf'};

xRaySourceCount      = 4;
backgroundPhotnRate  = 5.9e4;
timeBucket           = 1e5; % 1e5 sec
detectorArea         = 1; % m^2

earthEphemeris = loadEphemeris('earth', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
sunEphemeris   = loadEphemeris('sun', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
xRaySources    = loadXRaySources(xRaySourceCount);

initialXRay = loadInitialOrbit();
initialXRay = initialXRay(1:6);

% simulate real trajectory of spaceship
simulator          = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeDataXRay);
spaceshipTrueState = simulator.simulate(initialXRay);

xRayDetectorArgs.xRaySources = xRaySources;
xRayDetectorArgs.detectorArea = detectorArea;
xRayDetectorArgs.timeBucket = timeBucket;
xRayDetectorArgs.backgroundPhotnRate = backgroundPhotnRate;
xRayDetectorArgs.earthEphemeris = earthEphemeris;
xRayDetectorArgs.sunEphemeris = sunEphemeris;
xRayDetectorArgs.timeData = timeDataXRay;
xRayDetectorArgs.spaceshipState = spaceshipTrueState;


initArgsXRay.xRaySources = xRaySources;
initArgsXRay.earthEphemeris = [earthEphemeris.x(1); earthEphemeris.y(1); earthEphemeris.z(1)];
initArgsXRay.sunEphemeris = [sunEphemeris.x(1); sunEphemeris.y(1); sunEphemeris.z(1)];
initArgsXRay.invPeriods = getInvPeriods(xRaySources);
initArgsXRay.initialParams = [NaN NaN NaN];
initArgsXRay.observationNoiseMean = zeros(xRaySourceCount, 1);
initArgsXRay.observationNoiseCovariance = xRayToaCovariance(xRaySources, detectorArea, timeBucket, backgroundPhotnRate);
initArgsXRay.stateNoiseMean = [zeros(3, 1); zeros(3, 1)];
initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (3.25e-5*eye(3)).^2];

estimatorType = filterType(1);
iterations = zeros(iterationNumber, 2, timeDataXRay.SimulationNumber);

% parfor j = 1:iterationNumber
for j = 1:iterationNumber
    initialXRayCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (1.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    %initialXRayCov = [(30)^2*eye(3), zeros(3, 3); zeros(3, 3), (1e-5)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    initialXRayState = initialXRay + svdDecomposition(initialXRayCov)*randn(6, 1);
    
    trueState = spaceshipTrueState;
    xRayDetector  = XRayDetector(xRayDetectorArgs);
    xRayDetector.toa(iterationNumber == 1); % test, show X-Ray source signals
    
    xRayNavSystem = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, timeDataXRay, initArgsXRay, xRayDetector);
    
    stateEstimation   = xRayNavSystem.resolve(initialXRayState, 2*initialXRayCov, estimatorType, iterationNumber == 1);
    errTraj = vectNormError(trueState(1:3, :), stateEstimation(1:3, :), 1e3);
    errVel  = vectNormError(trueState(4:6, :), stateEstimation(4:6, :), 1e3);
    iterations(j, :, :) = [errTraj; errVel];
    
    if iterationNumber == 1 && length(estimatorType) == 1
        figure();
        subplot(2, 1, 1);
        e1 = (trueState(1:3, :) - stateEstimation(1:3, :))*1e3;
        plot2(timeDataXRay.Time, e1, 'trajectory error', {'x', 'y', 'z'}, 'trajectory error, m');
        subplot(2, 1, 2);
        e2 = (trueState(4:6, :) - stateEstimation(4:6, :))*1e3;
        plot2(timeDataXRay.Time, e2, 'velocity error', {'x', 'y', 'z'}, 'velocity error, m/sec');
    end
    
    fprintf('iteration of %d: completed\n', j );
end

if iterationNumber > 1
    errors = zeros(2, timeDataXRay.SimulationNumber);
    
    for i = 1:timeDataXRay.SimulationNumber
        errors(:, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    
    figure();
    subplot(2, 1, 1);
    plot2(timeDataXRay.Time, errors(1, :), 'trajectory errors', {'X-Ray Nav'}, 'trajectory error, meter');
    subplot(2, 1, 2);
    plot2(timeDataXRay.Time, errors(2, :), 'velocity errors', {'X-Ray Nav'}, 'velocity error, meter / sec');
end
