close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd   = '05:00:00.000';
timeData  = TimeExt(timeStart, timeEnd, 1e1, date, 1e5);

mass                        = 200;
iterationNumber             = 1;
secondInOneMinute           = 60;
esitimatedParams            = 2;
logLastErrors               = 1;
accelerationSigma           = 1e-7*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-7*ones(3, 1);      % [rad/sec]
accelerationInBodyFrame     = AccelerationInBodyFrame(timeData, [0; 0; 0], accelerationSigma);
angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeData, [0; 0; 0], angularVelocitySigma);
visualize                   = 1;

initialOrbit = loadInitialOrbit();

% debug free fly solver
simulator          = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeData, mass);
state = simulator.simulate(initialOrbit(1:6), visualize);

% debug controlled fly solver
simulator2 = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeData, mass);
stae2 = simulator2.simulate(visualize);
