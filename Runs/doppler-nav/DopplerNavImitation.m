close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd   = '00:00:50.000';
sunMoonInfluenceRefreshTime = 1; % sec
sampleTime         = 5e-3; % sec
timeDataDoppler    = TimeExt(timeStart, timeEnd, sampleTime, date, sunMoonInfluenceRefreshTime);
iterationNumber    = 120;
secondInOneMinute  = 120;
logData = 1;
initialOrbit = loadInitialOrbit();
initialDopplerNavState  = initialOrbit(1:6);

initialAcceleration     = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity  = zeros(3, 1);          % [rad/sec]
initialQuaternion       = initialOrbit(7:10);   % [-]
accelerationSigma       = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma    = 1e-4*ones(3, 1);      % [rad/sec]
dmuVariance             = (1*1e-5)^2; % [ (km / sec)^2 ] data from IET Radar Sonar Navig., 2011, Vol. 5, Iss. 9, pp. 1010-1017

% filterTypes = {'ukf', 'cdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'sghqf', 'pf'};
filterTypes = {'ckf'};

earthEphemeris = loadEphemeris('earth', timeDataDoppler.SimulationNumber, secondInOneMinute/timeDataDoppler.SampleTime);
sunEphemeris   = loadEphemeris('sun', timeDataDoppler.SimulationNumber, secondInOneMinute/timeDataDoppler.SampleTime);
accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataDoppler, initialAcceleration, accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataDoppler, initialAngularVelocity, angularVelocitySigma);

simulator          = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataDoppler);
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
initArgsDoppler.stateNoiseCovariance = [(2e0*eye(3)).^2 zeros(3); zeros(3) (2e-5*eye(3)).^2];

errors = zeros(length(filterTypes), 2, timeDataDoppler.SimulationNumber);
initialConditionStabilityKoeff = 1;
for l = 1:length(filterTypes)
    estimatorType = filterTypes(l);
    fprintf('estimator: %s\n', estimatorType{1});
    iterations = zeros(iterationNumber, 2, timeDataDoppler.SimulationNumber);
    estimationMatrix = zeros(iterationNumber, 6, timeDataDoppler.SimulationNumber);
    state = spaceshipTrueState.State;
    
    parfor j = 1:iterationNumber
%     for j = 1:iterationNumber
        initialDopplerCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (1.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
        initialDopplerState = initialDopplerNavState + initialConditionStabilityKoeff*chol(initialDopplerCov, 'lower')*randn(6, 1);
        
        trueState = state;
        k = 1e0; %1e2;
        
        dmu = DopplerMeasurementUnit(earthEphemeris, sunEphemeris, trueState, timeDataDoppler, dmuVariance);
        dopplerNavSystem = DopplerNavSystem(dmu, timeDataDoppler, initArgsDoppler, earthEphemeris, sunEphemeris);
        stateEstimation  = dopplerNavSystem.resolve(initialDopplerState, k*initialDopplerCov, estimatorType, iterationNumber == 1);
        
        errTraj = vectNormError(trueState(1:3, :), stateEstimation(1:3, :), 1);
        errVel  = vectNormError(trueState(4:6, :), stateEstimation(4:6, :), 1);
        iterations(j, :, :) = [errTraj; errVel];
        estimationMatrix(j, :, :) = [(trueState(1:3, :) - stateEstimation(1:3, :)); (trueState(1:3, :) - stateEstimation(1:3, :))];
                
        if iterationNumber == 1
            figure();
            subplot(2, 1, 1);
            e1 = (trueState(1:3, :) - stateEstimation(1:3, :))*1e3;
            plot2(timeDataDoppler.Time, e1, 'trajectory error', {'x', 'y', 'z'}, 'trajectory error, m');
            subplot(2, 1, 2);
            e2 = (trueState(4:6, :) - stateEstimation(4:6, :))*1e3;
            plot2(timeDataDoppler.Time, e2, 'velocity error', {'x', 'y', 'z'}, 'velocity error, m/sec');
        end
        
        fprintf('iteration of %d: completed\n', j );
    end
    
    % drow trajectory coord errors for each simulation
    figure();
    for i = 1:iterationNumber
        subplot(3, 1, 1);
        e1 = squeeze(estimationMatrix(i, 1, :));
        plot2(timeDataDoppler.Time, e1, 'trajectory error', {'r_x'}, 'trajectory error, km');
        
        subplot(3, 1, 2);
        e1 = squeeze(estimationMatrix(i, 2, :));
        plot2(timeDataDoppler.Time, e1, 'trajectory error', {'r_y'}, 'trajectory error, km');
        
        subplot(3, 1, 3);
        e1 = squeeze(estimationMatrix(i, 3, :));
        plot2(timeDataDoppler.Time, e1, 'trajectory error', {'r_z'}, 'trajectory error, km');
    end
    
    % drow velocity coord errors for each simulation
    figure();
    for i = 1:iterationNumber
        subplot(3, 1, 1);
        e1 = squeeze(estimationMatrix(i, 4, :));
        plot2(timeDataDoppler.Time, e1, 'velocity error', {'v_x'}, 'velocity error, km/sec');
        
        subplot(3, 1, 2);
        e1 = squeeze(estimationMatrix(i, 5, :));
        plot2(timeDataDoppler.Time, e1, 'velocity error', {'v_y'}, 'velocity error, km/sec');
        
        subplot(3, 1, 3);
        e1 = squeeze(estimationMatrix(i, 6, :));
        plot2(timeDataDoppler.Time, e1, 'velocity error', {'v_z'}, 'velocity error, km/sec');
    end
    
    if iterationNumber > 1
        for i = 1:timeDataDoppler.SimulationNumber
            errors(l, :, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
        end
    end
end

if iterationNumber > 1
    %     fileName = strcat('errors_dns.mat');
    %     save(fileName, 'errors');
    
    figure();
    subplot(2, 1, 1);
    plot2(timeDataDoppler.Time, squeeze(errors(:, 1, :)), 'trajectory errors', filterTypes, 'trajectory error, km');
    subplot(2, 1, 2);
    plot2(timeDataDoppler.Time, squeeze(errors(:, 2, :)), 'velocity errors', filterTypes, 'velocity error, km / sec');
    
    for i = 1:length(filterTypes)
        estimatorType = filterTypes(i);
        fprintf('estimator: %s\n', estimatorType{1});
        fprintf('trajectory RMS: %d\n', squeeze(errors(i, 1, end)));
        fprintf('velocity RMS: %d\n', squeeze(errors(i, 2, end)));
    end
end
