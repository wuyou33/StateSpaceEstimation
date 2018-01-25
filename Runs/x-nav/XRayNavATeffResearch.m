close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);

date.day            = 17;
date.mon            = 11;
date.year           = 2017;
timeStart           = '00:00:00.000';
timeEnd             = '01:30:00.000';
timeDataXRay        = TimeExt(timeStart, timeEnd, 100, date, 1e7); % change refreshSunMoonInfluenceTime to real number
iterationNumber     = 700;
secondInOneMinute   = 60;
esitimatedParams    = 2;
logLastErrors       = 1;
mass                = 200; % [kg]
errorBudget         = 20; % [%]
timeBucketArray     = [5e2 1e3 5e3 1e4 5e4 1e5];
timeBucketLegend    = {'AT_e_f_f = 5e2', 'AT_e_f_f = 1e3', 'AT_e_f_f = 5e3', 'AT_e_f_f = 1e4', 'AT_e_f_f = 5e4', 'AT_e_f_f = 1e5'};
%{'ukf', 'srukf', 'cdkf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'ekf', 'pf', 'sppf', 'gspf', 'gmsppf'};
estimatorType = {'srukf'};

set(0, 'defaultfigurecolor', [1 1 1]);

b_det   = 0.1; % Detector Background Rate. [photon*cm^2*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^2*sec^-1]
xRaySourceCount      = 4;
backgroundPhotonRate = b_det + b_diff + b_cosm;
detectorArea         = 1; % [m^2]

earthEphemeris = loadEphemeris('earth', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
sunEphemeris   = loadEphemeris('sun', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
xRaySources    = loadXRaySources(xRaySourceCount);

initialXRay = loadInitialOrbit();
initialXRay = initialXRay(1:6);

% simulate real trajectory of spaceship
simulator = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeDataXRay, mass);
trueState = simulator.simulate(initialXRay, iterationNumber == 1);

errors = zeros(length(timeBucketArray), esitimatedParams, timeDataXRay.SimulationNumber);
e1 = zeros(3, timeDataXRay.SimulationNumber);
e2 = zeros(3, timeDataXRay.SimulationNumber);

errLegend = {1:length(timeBucketArray)};
tic
for l = 1:length(timeBucketArray)
    timeBucket = timeBucketArray(l);
    
    xRayDetectorArgs.xRaySources = xRaySources;
    xRayDetectorArgs.detectorArea = detectorArea;
    xRayDetectorArgs.timeBucket = timeBucket;
    xRayDetectorArgs.backgroundPhotonRate = backgroundPhotonRate;
    xRayDetectorArgs.earthEphemeris = earthEphemeris;
    xRayDetectorArgs.sunEphemeris = sunEphemeris;
    xRayDetectorArgs.timeData = timeDataXRay;
    xRayDetectorArgs.spaceshipState = trueState;
    xRayDetectorArgs.errorBudget = errorBudget;
    
    
    initArgsXRay.xRaySources = xRaySources;
    initArgsXRay.earthEphemeris = [earthEphemeris.x(1); earthEphemeris.y(1); earthEphemeris.z(1)];
    initArgsXRay.sunEphemeris = [sunEphemeris.x(1); sunEphemeris.y(1); sunEphemeris.z(1)];
    initArgsXRay.invPeriods = getInvPeriods(xRaySources);
    initArgsXRay.initialParams = [NaN NaN NaN];
    initArgsXRay.observationNoiseMean = zeros(xRaySourceCount, 1);
    initArgsXRay.observationNoiseCovariance = xRayToaCovariance(xRaySources, detectorArea, timeBucket, backgroundPhotonRate, errorBudget);
    initArgsXRay.stateNoiseMean = [zeros(3, 1); zeros(3, 1)];
    initArgsXRay.mass = mass;
    initArgsXRay.startTime = timeDataXRay.StartSecond;
    initArgsXRay.gravityModel = m_fitSolarSystemGravityModel(timeDataXRay.SampleTime, timeDataXRay.SimulationNumber);
    
    initialXRayCov = [(5)^2*eye(3), zeros(3, 3); zeros(3, 3), (5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    %     initialXRayCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (0.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    
    if stringmatch(estimatorType, {'ekf'})
        initArgsXRay.stateNoiseCovariance = [(9.5e-1*eye(3)).^2 zeros(3); zeros(3) (5e-2*eye(3)).^2];
    elseif stringmatch(estimatorType, {'gmsppf'})
        initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (1e-4*eye(3)).^2];
    elseif stringmatch(estimatorType, {'sppf'})
        initArgsXRay.stateNoiseCovariance = [(7.5e-4*eye(3)).^2 zeros(3); zeros(3) (1e-4*eye(3)).^2];
    else
        initArgsXRay.stateNoiseCovariance = [(1e-2*eye(3)).^2 zeros(3); zeros(3) (8.5e-4*eye(3)).^2];
    end
    
    alpha   = (1)^2;
    beta    = (1)^2;
    gamma   = (1)^2;
    
    initArgsXRay.stateNoiseCovariance = beta * initArgsXRay.stateNoiseCovariance;
    initArgsXRay.observationNoiseCovariance = gamma * initArgsXRay.observationNoiseCovariance;
    
    iterations = zeros(iterationNumber, 2, timeDataXRay.SimulationNumber);
    x_iterations = zeros(iterationNumber, 6, timeDataXRay.SimulationNumber);
    
    parfor j = 1:iterationNumber
        tstate = trueState;
        initialXRayState = initialXRay + chol(initialXRayCov, 'lower') * randn(6, 1);
        xRayDetector  = XRayDetector(xRayDetectorArgs);
        xRayDetector.toa(iterationNumber == 1); % test, show X-Ray source signals
        
        xRayNavSystem = XRayNavSystem(earthEphemeris, sunEphemeris, xRaySources, timeDataXRay, initArgsXRay, xRayDetector);
        stateEstimation = xRayNavSystem.resolve(initialXRayState, alpha*initialXRayCov, estimatorType, iterationNumber == 1);
        
        errTraj = vectNormError(tstate(1:3, :), stateEstimation(1:3, :), 1e3);
        errVel  = vectNormError(tstate(4:6, :), stateEstimation(4:6, :), 1e3);
        
        iterations(j, :, :) = [errTraj; errVel];
        x_iterations(j, :, :) = stateEstimation;
    end
    
    for i = 1:timeDataXRay.SimulationNumber
        errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    
    fprintf('timeBucket: %d\n', timeBucketArray(l));
    fprintf('RMS trajectory: %d\n', errors(l, 1, end));
    fprintf('RMS velocity: %d\n', errors(l, 2, end));
    
    fprintf('imitation # %d from %d completed\n', l,  length(timeBucketArray));
end
toc

figure();
subplot(2, 1, 1);
plot2(timeDataXRay.TimeInHour, squeeze(errors(:, 1, :)), 'RMS trajectory error (X-Ray nav)', timeBucketLegend, 'RMS trajectory, meter', 'time, hours');
hold on;
subplot(2, 1, 2);
plot2(timeDataXRay.TimeInHour, squeeze(errors(:, 2, :)), 'RMS velocity error (X-Ray nav)', timeBucketLegend, 'RMS velocity, meter / sec', 'time, hours');
hold on;