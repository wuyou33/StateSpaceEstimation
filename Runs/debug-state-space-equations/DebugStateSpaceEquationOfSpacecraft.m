close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

date.day  = 17;
date.mon  = 11;
date.year = 2015;
timeStart = '00:00:00.000';
timeEnd   = '05:00:00.000';
timeData  = TimeExt(timeStart, timeEnd, 1, date, 1e5);

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
state2 = simulator2.simulate(visualize);


insInitArgs.accBiasMu                    = zeros(3, 1);
insInitArgs.accBiasSigma                 = zeros(3, 1);
insInitArgs.gyroBiasMu                   = zeros(3, 1);
insInitArgs.gyroBiasSigma                = zeros(3, 1);
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.visualize                    = 1;
insInitArgs.timeData                     = timeData;
insInitArgs.accNoiseVar                  = zeros(3, 1);
insInitArgs.gyroNoiseVar                 = zeros(3, 1);
insInitArgs.accScale                     = zeros(3);
insInitArgs.gyroScale                    = zeros(3);
insInitArgs.levelArm                     = zeros(3, 1);
insInitArgs.angularAccelerBodyFrame      = zeros(3, 1);
insInitArgs.gyroGSensitiveBias           = zeros(3);
insInitArgs.gravityModel                 = fitSolarSystemGravityModel(timeData.SampleTime, timeData.SimulationNumber);
insInitArgs.mass                         = mass;

ins = initInertialNavigationSystem('init', insInitArgs);
estimations = ins.simulate(initialOrbit);

figure();
r_err = (state2.Trajectory - estimations(1:3, :))';
plot(r_err); grid on;

figure();
v_err = (state2.Velocity - estimations(4:6, :))';
plot(v_err); grid on;

figure();
q_err = (state2.FullState(7:10, :) - estimations(7:10, :))';
plot(q_err); grid on;
