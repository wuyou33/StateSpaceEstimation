close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd = '00:00:00.500';
sunMoonInfluenceRefreshTime = 1; % sec
sampleTime = 1e-3; % sec
timeDataXRayDoppler  = TimeExt(timeStart, timeEnd, sampleTime, date, sunMoonInfluenceRefreshTime);

iterationNumber    = 1;
secondInOneMinute  = 60;

initialOrbit = load_initial_orbit();
initialXRayDopplerNavState = initialOrbit(1:6);

initialAcceleration     = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity  = zeros(3, 1);          % [rad/sec]
initialQuaternion       = initialOrbit(7:10);   % [-]
accelerationSigma       = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma    = 1e-4*ones(3, 1);      % [rad/sec]
dmuVariance             = (1e-5)^2; % [ (km / sec)^2 ] data from IET Radar Sonar Navig., 2011, Vol. 5, Iss. 9, pp. 1010-1017
xRaySourceCount         = 4;
backgroundPhotnRate     = 5.9e4;
timeBucket              = 1e5; % 1e5 sec
detectorArea            = 1; % m^2

xRaySources             = load_x_ray_sources(xRaySourceCount);

%{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf', 'sppf', 'fdckf', 'fdckfAugmented', 'cqkf', 'gspf', 'gmsppf', 'ghqf', 'sghqf'};
filterTypeXRay    = {'ckf'};
filterTypeDoppler = {'pf'};

earthEphemeris = load_ephemeris('earth', timeDataXRayDoppler, secondInOneMinute/timeDataXRayDoppler.SampleTime);
sunEphemeris   = load_ephemeris('sun', timeDataXRayDoppler, secondInOneMinute/timeDataXRayDoppler.SampleTime);
accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataXRayDoppler, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataXRayDoppler, initialAngularVelocity, angularVelocitySigma);

simulator          = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataXRayDoppler);
spaceshipTrueState = simulator.simulate(iterationNumber == 1);

initArgsDoppler.earthEphemeris.x = earthEphemeris.x(1);
initArgsDoppler.earthEphemeris.y = earthEphemeris.y(1);
initArgsDoppler.earthEphemeris.z = earthEphemeris.z(1);
initArgsDoppler.earthEphemeris.vx = earthEphemeris.vx(1);
initArgsDoppler.earthEphemeris.vy = earthEphemeris.vy(1);
initArgsDoppler.earthEphemeris.vz = earthEphemeris.vz(1);

initArgsDoppler.sunEphemeris.x = sunEphemeris.x(1);
initArgsDoppler.sunEphemeris.y = sunEphemeris.y(1);
initArgsDoppler.sunEphemeris.z = sunEphemeris.z(1);
initArgsDoppler.sunEphemeris.vx = sunEphemeris.vx(1);
initArgsDoppler.sunEphemeris.vy = sunEphemeris.vy(1);
initArgsDoppler.sunEphemeris.vz = sunEphemeris.vz(1);

initArgsDoppler.initialParams = [NaN NaN NaN];

initArgsDoppler.observationNoiseMean = 0;
initArgsDoppler.observationNoiseCovariance = dmuVariance;
initArgsDoppler.stateNoiseMean = [zeros(3, 1); zeros(3, 1)];
initArgsDoppler.stateNoiseCovariance = [(2e-5*eye(3)).^2 zeros(3); zeros(3) (2e-7*eye(3)).^2];


xRayDetectorArgs.xRaySources = xRaySources;
xRayDetectorArgs.detectorArea = detectorArea;
xRayDetectorArgs.timeBucket = timeBucket;
xRayDetectorArgs.backgroundPhotnRate = backgroundPhotnRate;
xRayDetectorArgs.earthEphemeris = earthEphemeris;
xRayDetectorArgs.sunEphemeris = sunEphemeris;
xRayDetectorArgs.timeData = timeDataXRayDoppler;
xRayDetectorArgs.spaceshipState = spaceshipTrueState.State;


initArgsXRay.xRaySources = xRaySources;
initArgsXRay.earthEphemeris = [earthEphemeris.x(1); earthEphemeris.y(1); earthEphemeris.z(1)];
initArgsXRay.sunEphemeris = [sunEphemeris.x(1); sunEphemeris.y(1); sunEphemeris.z(1)];
initArgsXRay.invPeriods = get_inv_periods(xRaySources);
initArgsXRay.initialParams = [NaN NaN NaN];
initArgsXRay.observationNoiseMean = zeros(xRaySourceCount, 1);
initArgsXRay.observationNoiseCovariance = xRayToaCovariance(xRaySources, detectorArea, timeBucket, backgroundPhotnRate);
initArgsXRay.stateNoiseMean = [zeros(3, 1); zeros(3, 1)];
initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (3.25e-5*eye(3)).^2];


state = spaceshipTrueState.State;
iterations = zeros(iterationNumber, 2, timeDataXRayDoppler.SimulationNumber);
report = iterationNumber == 1;

% parfor j = 1:iterationNumber
for j = 1:iterationNumber
    initXRayCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (1.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    %initialXRayCov = [(30)^2*eye(3), zeros(3, 3); zeros(3, 3), (1e-5)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    initXRayState = initialXRayDopplerNavState + chol(initXRayCov, 'lower')*randn(6, 1);
    
    initDopplerCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (1.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    initDopplerState = initialXRayDopplerNavState + chol(initDopplerCov, 'lower')*randn(6, 1);
    
    trueState = state;
    
    dmu = DopplerMeasurementUnit(earthEphemeris, sunEphemeris, trueState, timeDataXRayDoppler, dmuVariance);
    dopplerNavSystem = DopplerNavSystem(dmu, timeDataXRayDoppler, initArgsDoppler, earthEphemeris, sunEphemeris);
    
    xRayDetector  = XRayDetector(xRayDetectorArgs);        
    x_ray_nav_system = X_RayNavigationSystem(earthEphemeris, sunEphemeris, xRaySources, timeDataXRayDoppler, initArgsXRay, xRayDetector);    

    xRayDoppNS = XRayDopplerIntegrated(timeDataXRayDoppler, x_ray_nav_system, dopplerNavSystem);    
    estimations = xRayDoppNS.resolve(filterTypeXRay(1), filterTypeDoppler(1), initXRayState, initXRayCov, initDopplerState, initDopplerCov, report);
    
    errTraj = vect_norm_error(trueState(1:3, :), estimations(1:3, :), 1e3);
    errVel  = vect_norm_error(trueState(4:6, :), estimations(4:6, :), 1e3);
    iterations(j, :, :) = [errTraj; errVel];
    
    if iterationNumber == 1
        figure();
        subplot(2, 1, 1);
        e1 = (trueState(1:3, :) - estimations(1:3, :))*1e3;
        plot2(timeDataXRayDoppler.Time, e1, 'trajectory error', {'x', 'y', 'z'}, 'trajectory error, m');
        subplot(2, 1, 2);
        e2 = (trueState(4:6, :) - estimations(4:6, :))*1e3;
        plot2(timeDataXRayDoppler.Time, e2, 'velocity error', {'x', 'y', 'z'}, 'velocity error, m/sec');
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
    plot2(timeDataXRayDoppler.Time, errors(1, :), 'trajectory errors', {'X-Ray Nav'}, 'trajectory error, meter');
    subplot(2, 1, 2);
    plot2(timeDataXRayDoppler.Time, errors(2, :), 'velocity errors', {'X-Ray Nav'}, 'velocity error, meter / sec');
end
