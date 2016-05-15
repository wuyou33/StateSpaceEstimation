close all force; clc;

% 
% anonymousFunc = @(x, y) [...
%     x^2;
%     x*y...
%     ];

addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('InertialMeasurementUnit')); % include folder with inertial measurement functions
addpath(genpath('InertialNavigationSystem')); % include folder with inertial navigation system functions
addpath(genpath('SNS')); % include folder with inertial navigation system functions
addpath(genpath('Integrated')); % include folder with integrated navigation system functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('TOAMeasurementUnit')); % include folder with Pulsar & Quasar navigation
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));
addpath(genpath('XNAV'));

filterTypeArray         = {'sckf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf'};

sampleTime              = 1; % seconds
simulationTime          = 1*20*60; % y hours * 60 min * 60 sec
simulationNumber        = round(simulationTime / sampleTime);
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;

accBiasMu               = 1e-8*[1;1;1];
accBiasSigma            = 5e-7*[1;1;1];
accNoiseVar             = 1e-8*[1;1;1];

gyroBiasMu              = 1e-8*[1;1;1];
gyroBiasSigma           = 5e-6*[1;1;1];
gyroNoiseVar            = 1e-5*[1;1;1];

initialAcceleration     = zeros(3,1);
initialAngularVelocity  = zeros(3,1);
initialQuaternion       = [1;0;0;0];

%%  
accelerationInBodyFrame = AccelerationInBodyFrame(simulationNumber, ... simulationNumber
    sampleTime, ... sampleTime
    [0; 0; 0], ... mu
    1.75*1e-5*[1;1;1] ... sigma    
);

angularVelocityInBodyFrame = AngularVelocityInBodyFrame(simulationNumber, ... simulationNumber
    sampleTime, ... sampleTime
    [0; 0; 0], ... mu
    1.2*1e-3*[1;1;1] ... sigma
);

T_till_current_epoch = 0.1465;
initialSpaceshipState = loadInitialOrbit();

tic;

simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialSpaceshipState, ...
    accelerationInBodyFrame, ... accelerationInBodyFrame
    angularVelocityInBodyFrame, ... angularVelocityInBodyFrame
    T_till_current_epoch);

spaceshipStateTrue = simulator.Simulate(time, sampleTime);

fprintf('Simulating true trajectory: ');
toc;

SatelliteOrbitVisualization(spaceshipStateTrue);

%% simulate INS
satellitePhaseState = SatellitePhaseSpace(initialSpaceshipState, simulationNumber);
iterations = repmat(satellitePhaseState, iterationNumber, 1);

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.simulationNumber             = simulationNumber;
insInitArgs.timeMinutes                  = timeMinutes;
insInitArgs.visualize                    = 0; % 1
insInitArgs.T_till_current_epoch         = T_till_current_epoch;
insInitArgs.sampleTime                   = sampleTime;
insInitArgs.accNoiseVar                  = accNoiseVar;
insInitArgs.gyroNoiseVar                 = gyroNoiseVar;

snsInitArgs.Trajectory = spaceshipStateTrue.Trajectory';
snsInitArgs.Velocity   = spaceshipStateTrue.Velocity';

for j = 1:length(filterTypeArray)
    estimatorType = filterTypeArray(j);
    fprintf('estimator: %s\n', estimatorType{1});
    tic;
    
%     parfor k = 1:iterationNumber    
    for k = 1:iterationNumber    
        insSns = IntegratedInsSns(insInitArgs, snsInitArgs, initialAcceleration, initialAngularVelocity, initialQuaternion, time, sampleTime);        
        iterations(k) = insSns.simulate(initialSpaceshipState + [1500*randn(3,1); 2.75*randn(3,1); 0.0001*randn(4,1)], estimatorType, 1, spaceshipStateTrue);
        
        fprintf ('thread of %d: ', k );
    end

    fprintf('Simulation: ');
    toc;    
end