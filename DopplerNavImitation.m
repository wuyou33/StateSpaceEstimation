close all; clc; clearvars;

addpath(genpath('./'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd = '00:00:02.000';
sunMoonInfluenceRefreshTime = 1; % sec
sampleTime = 1e-1; % sec
timeDataDoppler  = TimeExt(timeStart, timeEnd, sampleTime, date, sunMoonInfluenceRefreshTime);

iterationNumber    = 1;
secondInOneMinute  = 60;

initialOrbit = loadInitialOrbit();
initialDopplerNavState  = initialOrbit(1:6);

initialAcceleration     = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity  = zeros(3, 1);          % [rad/sec]
initialQuaternion       = initialOrbit(7:10);   % [-]
accelerationSigma       = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma    = 1e-4*ones(3, 1);      % [rad/sec]
dmuVariance             = (1e-5)^2; % [ (km / sec)^2 ] data from IET Radar Sonar Navig., 2011, Vol. 5, Iss. 9, pp. 1010-1017

filterType = {'ukf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf', 'sppf', 'fdckf', 'fdckfAugmented', 'cqkf', 'gspf', 'gmsppf', 'ghqf', 'sghqf'};

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
initArgsDoppler.stateNoiseCovariance = [(2e-5*eye(3)).^2 zeros(3); zeros(3) (2e-7*eye(3)).^2];

estimatorType = filterType(1);
iterations = zeros(iterationNumber, 2, timeDataDoppler.SimulationNumber);

state = spaceshipTrueState.State;
% parfor j = 1:iterationNumber
for j = 1:iterationNumber
    initialDopplerCov = [(0.5)^2*eye(3), zeros(3, 3); zeros(3, 3), (1.5e-3)^2*eye(3)]; % [ [km^2], [(km / sec)^2] ]
    initialDopplerState = initialDopplerNavState + svdDecomposition(initialDopplerCov)*randn(6, 1);
    
    trueState = state;
    
    dmu = DopplerMeasurementUnit(earthEphemeris, sunEphemeris, trueState, timeDataDoppler, dmuVariance);
    % dmu.simulate(iterationNumber == 1);
    
    tic;
    dopplerNavSystem = DopplerNavSystem(dmu, timeDataDoppler, initArgsDoppler, earthEphemeris, sunEphemeris);
    stateEstimation   = dopplerNavSystem.resolve(initialDopplerState, 2*initialDopplerCov, estimatorType, iterationNumber == 1);
    toc;
    
    errTraj = vectNormError(trueState(1:3, :), stateEstimation(1:3, :), 1);
    errVel  = vectNormError(trueState(4:6, :), stateEstimation(4:6, :), 1);
    iterations(j, :, :) = [errTraj; errVel];
    
    if iterationNumber == 1 && length(estimatorType) == 1
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

if iterationNumber > 1
    errors = zeros(2, timeDataDoppler.SimulationNumber);
    
    for i = 1:timeDataDoppler.SimulationNumber
        errors(:, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
    end
    
    figure();
    subplot(2, 1, 1);
    plot2(timeDataDoppler.Time, errors(1, :), 'trajectory errors', {'Doppler Nav'}, 'trajectory error, km');
    subplot(2, 1, 2);
    plot2(timeDataDoppler.Time, errors(2, :), 'velocity errors', {'Doppler Nav'}, 'velocity error, km / sec');
end
