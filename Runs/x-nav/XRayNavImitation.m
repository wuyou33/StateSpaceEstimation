close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

% RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));

m_fit_solar_system_gravity_model = memoize(@fit_solar_system_gravity_model);

date.day            = 17;
date.mon            = 11;
date.year           = 2017;
timeStart           = '00:00:00.000';
timeEnd             = '05:00:00.000';
timeDataXRay        = TimeExt(timeStart, timeEnd, 100, date, 1e7); % change refreshSunMoonInfluenceTime to real number
iterationNumber     = 1;
secondInOneMinute   = 60;
esitimatedParams    = 2;
drawIterations      = 1;
mass                = 200; % [kg]
errorBudget         = 20; % [%]

%{'ukf', 'srukf', 'cdkf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'ekf', 'pf', 'sppf', 'gspf', 'gmsppf'};
filterTypes = {'srukf'};
 

b_det   = 0.1; % Detector Background Rate. [photon*cm^2*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^2*sec^-1]
xRaySourceCount      = 4;
backgroundPhotonRate = b_det + b_diff + b_cosm;
timeBucket           = 1e5; % [sec]
detectorArea         = 1; % [m^2]

earthEphemeris = load_ephemeris('earth', timeDataXRay, secondInOneMinute/timeDataXRay.SampleTime);
sunEphemeris   = load_ephemeris('sun', timeDataXRay, secondInOneMinute/timeDataXRay.SampleTime);
xRaySources    = load_x_ray_sources(xRaySourceCount);

initialXRay = load_initial_orbit();
initialXRay = initialXRay(1:6);

% simulate real trajectory of spaceship
simulator = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeDataXRay, mass);
trueState = simulator.simulate(initialXRay, iterationNumber == 1);

errors = zeros(length(filterTypes), esitimatedParams, timeDataXRay.SimulationNumber);
e1 = zeros(3, timeDataXRay.SimulationNumber);
e2 = zeros(3, timeDataXRay.SimulationNumber);

for l = 1:length(filterTypes)
    estimatorType = filterTypes(l);
    
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
    initArgsXRay.invPeriods = get_inv_periods(xRaySources);
    initArgsXRay.initialParams = [NaN NaN NaN];
    initArgsXRay.observationNoiseMean = zeros(xRaySourceCount, 1);
    initArgsXRay.observationNoiseCovariance = xRayToaCovariance(xRaySources, detectorArea, timeBucket, backgroundPhotonRate, errorBudget);
    initArgsXRay.stateNoiseMean = [zeros(3, 1); zeros(3, 1)];
    initArgsXRay.mass = mass;
    initArgsXRay.startTime = timeDataXRay.StartSecond;
    initArgsXRay.gravityModel = m_fit_solar_system_gravity_model(timeDataXRay.SampleTime, timeDataXRay.SimulationNumber);
    
	initialXRayCov = [(5)^2*eye(3), zeros(3, 3); zeros(3, 3), (5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
%     initialXRayCov = [(5)^2*eye(3), zeros(3, 3); zeros(3, 3), (0.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    
    if string_match(estimatorType, {'pf', 'gspf'})
        initArgsXRay.stateNoiseCovariance = [(1e-2*eye(3)).^2 zeros(3); zeros(3) (5e-5*eye(3)).^2];
    elseif string_match(estimatorType, {'gmsppf'})
        initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (1e-4*eye(3)).^2];
    elseif string_match(estimatorType, {'sppf'})
        initArgsXRay.stateNoiseCovariance = [(1e-3*eye(3)).^2 zeros(3); zeros(3) (8.5e-5*eye(3)).^2];
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
    tic;
    fprintf('estimator: %s\n', estimatorType{1});
    for j = 1:iterationNumber
        rng(sum(100*clock), 'v5normal');
        rng(sum(100*clock), 'v5uniform');
        try
            tstate = trueState;
            initialXRayState = initialXRay + chol(initialXRayCov, 'lower') * randn(6, 1);
            xRayDetector  = XRayDetector(xRayDetectorArgs);
            xRayDetector.toa(iterationNumber == 1); % test, show X-Ray source signals
            
            x_ray_nav_system = X_RayNavigationSystem(earthEphemeris, sunEphemeris, xRaySources, timeDataXRay, initArgsXRay, xRayDetector);
            stateEstimation = x_ray_nav_system.resolve(initialXRayState, alpha*initialXRayCov, estimatorType, iterationNumber == 1);
            
            errTraj = vect_norm_error(tstate(1:3, :), stateEstimation(1:3, :), 1e3);
            errVel  = vect_norm_error(tstate(4:6, :), stateEstimation(4:6, :), 1e3);
            
            iterations(j, :, :) = [errTraj; errVel];
            x_iterations(j, :, :) = stateEstimation;
            
            fprintf('iteration of %d: completed\n', j );
        catch exc
            display(exc.message);
            rethrow(exc)
        end
    end
    toc;
    
    if drawIterations
        figure();
        for j = 1:iterationNumber
            x_est = squeeze(x_iterations(j, :, :));
            
            subplot(2, 1, 1);
            e1 = (trueState(1:3, :) - x_est(1:3, :))*1e3;
            plot2(timeDataXRay.TimeInHour, e1, sprintf('Ensemble of trajectory estimation errors (%s)', estimatorType{1}), ...
                {'Trajectory error'}, 'Trajectory error, meter', 'time, hours');
            
            hold on;
            subplot(2, 1, 2);
            e2 = (trueState(4:6, :) - x_est(4:6, :))*1e3;
            plot2(timeDataXRay.TimeInHour, e2, sprintf('Ensemble of velocity estimation errors (%s)', estimatorType{1}), ...
                {'Velocity error'}, 'Velocity error, meter / sec', 'time, hours');
            hold on;
        end
    end
    
    if iterationNumber > 1
        for i = 1:timeDataXRay.SimulationNumber
            errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
        end
    end
    
    if drawIterations
        figure();
        subplot(2, 1, 1);
        e1 = (trueState(1:3, :) - x_est(1:3, :))*1e3;
        plot2(timeDataXRay.TimeInHour, e1, sprintf('Transient response of estimation of trajectory error (%s)', estimatorType{1}),...
            {'x', 'y', 'z'}, 'Trajectory error, meter', 'time, hours');
        subplot(2, 1, 2);
        e2 = (trueState(4:6, :) - x_est(4:6, :))*1e3;
        plot2(timeDataXRay.TimeInHour, e2, sprintf('Transient response of estimation of velocity error (%s)', estimatorType{1}), ...
            {'x', 'y', 'z'}, 'Velocity error, meter / sec', 'time, hours');
    end
end

if iterationNumber > 1
    %     save('errors_xray_sp_family.mat', 'errors');
    figure();
    subplot(2, 1, 1);
    plot2(timeDataXRay.TimeInHour, squeeze(errors(:, 1, :)), 'RMS trajectory error (X-Ray nav)', filterTypes, 'RMS trajectory, meter', 'time, hours');
    subplot(2, 1, 2);
    plot2(timeDataXRay.TimeInHour, squeeze(errors(:, 2, :)), 'RMS velocity error (X-Ray nav)', filterTypes, 'RMS velocity, meter / sec', 'time, hours');
    
    for i = 1:length(filterTypes)
        estimatorType = filterTypes(i);
        fprintf('estimator: %s\n', estimatorType{1});
        fprintf('RMS trajectory: %d\n', errors(i, 1, end));
        fprintf('RMS velocity: %d\n', errors(i, 2, end));
    end
end
