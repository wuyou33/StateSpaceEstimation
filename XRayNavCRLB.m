close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);

date.day            = 17;
date.mon            = 11;
date.year           = 2017;
timeStart           = '00:00:00.000';
timeEnd             = '00:20:00.000';
timeDataXRay        = TimeExt(timeStart, timeEnd, 1e1, date, 1e7); % change refreshSunMoonInfluenceTime to real number
iterationNumber     = 4;
secondInOneMinute   = 60;
esitimatedParams    = 2;
logLastErrors       = 1;
mass                = 200; % [kg]
errorBudget         = 10; % [%]

b_det   = 0.1; % Detector Background Rate. [photon*cm^2*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^2*sec^-1]
xRaySourceCount      = 4;
backgroundPhotonRate = b_det + b_diff + b_cosm;
timeBucket           = 1e5; % 1e5 sec
detectorArea         = 1; % m^2

earthEphemeris = loadEphemeris('earth', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
sunEphemeris   = loadEphemeris('sun', timeDataXRay.SimulationNumber, secondInOneMinute/timeDataXRay.SampleTime);
xRaySources    = loadXRaySources(xRaySourceCount);

initialXRay = loadInitialOrbit();
initialXRay = initialXRay(1:6);

% simulate real trajectory of spaceship
simulator          = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeDataXRay, mass);
trueState = simulator.simulate(initialXRay, iterationNumber == 1);

errors = zeros(esitimatedParams, timeDataXRay.SimulationNumber);
initialConditionStabilityKoeff = 1;


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
initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (1e-4*eye(3)).^2];
% initArgsXRay.stateNoiseCovariance = [(9.5e-2*eye(3)).^2 zeros(3); zeros(3) (5e-4*eye(3)).^2];

iterations = zeros(iterationNumber, 2, timeDataXRay.SimulationNumber - 1);

args.type  = 'state';
args.tag   = 'State estimation for X-Ray navigation system';
args.model = gssmXNav('init', initArgsXRay);
[stateNoise, observNoise, inferenceModel] = inferenceNoiseGenerator(inferenceDataGenerator(args), 'ekf');
initialXRayCov = [(5)^2*eye(3), zeros(3, 3); zeros(3, 3), (5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]

bound = zeros(6, timeDataXRay.SimulationNumber - 1);

tEpoch = currentEpoch(timeDataXRay.JD, timeDataXRay.StartSecond);
for j = 1:iterationNumber
    initialXRayState = initialXRay + initialConditionStabilityKoeff*chol(initialXRayCov, 'lower') * randn(6, 1);
    xRayDetector  = XRayDetector(xRayDetectorArgs);
    dtoa = xRayDetector.toa(0); % test, show X-Ray source signals
    crlb_predict = initialXRayCov;
    
    for i = 2:timeDataXRay.SimulationNumber
        modelParams(1) = tEpoch;
        modelParams(2) = timeDataXRay.SampleTime;
        modelParams(3) = timeDataXRay.Time(i);
        
        earthEphemerisStep = [earthEphemeris.x(i); earthEphemeris.y(i); earthEphemeris.z(i)];
        sunEphemerisStep   = [sunEphemeris.x(i); sunEphemeris.y(i); sunEphemeris.z(i)];
        
        updModel = inferenceModel.model.setParams(inferenceModel.model, ...
            modelParams, ...
            inferenceModel.model.xRaySources, ...
            earthEphemerisStep, ...
            sunEphemerisStep, ...
            inferenceModel.model.invPeriods, ...
            inferenceModel.model.mass, ...
            inferenceModel.model.gravityModel, ...
            inferenceModel.model.startTime);
        inferenceModel.model = updModel;
        [bound(:, i - 1), crlb_predict] = crlb(crlb_predict, trueState(:, i), dtoa(:, i), stateNoise, observNoise, inferenceModel);
    end
    
    errTraj = sum( (bound(1:3, :) * 1e3).^2, 1 ).^0.5;
    errVel  = sum( (bound(4:6, :) * 1e3).^2, 1 ).^0.5;
    iterations(j, :, :) = [errTraj; errVel];
    
    if iterationNumber == 1
        figure();
        subplot(2, 1, 1);
        plot2(timeDataXRay.TimeInHour(2:end), bound(1:3, :) * 1e3, 'Position error', {'x', 'y', 'z'}, 'CRLB position error, m', 'time, hours');
        subplot(2, 1, 2);
        plot2(timeDataXRay.TimeInHour(2:end), bound(4:6, :) * 1e3, 'Velocity error', {'x', 'y', 'z'}, 'CRLB velocity error, m/sec', 'time, hours');
    end
    
    fprintf('iteration of %d: completed\n', j );
end

if iterationNumber > 1
    for i = 1:timeDataXRay.SimulationNumber-1
        errors(:, i+1) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    errors(1, 1) = 1e3*norm([initialXRayCov(1, 1)^.5 initialXRayCov(2, 2)^.5 initialXRayCov(3, 3)^.5]);
    errors(2, 1) = 1e3*norm([initialXRayCov(4, 4)^.5 initialXRayCov(5, 5)^.5 initialXRayCov(6, 6)^.5]);
end

if iterationNumber > 1
    figure();
    subplot(2, 1, 1);
    plot2(timeDataXRay.TimeInHour, squeeze(errors(1, :)), 'CRLB position error in X-Ray', {'CRLB_r'}, 'CRLB trajectory error, meter', 'time, hours');
    subplot(2, 1, 2);
    plot2(timeDataXRay.TimeInHour, squeeze(errors(2, :)), 'CRLB velocity error in X-Ray', {'CRLB_v'}, 'CRLB velocity error, meter / sec', 'time, hours');
end