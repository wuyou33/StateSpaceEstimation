close all force; clc;

addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('InertialMeasurementUnit')); % include folder with inertial measurement functions
addpath(genpath('InertialNavigationSystem')); % include folder with inertial navigation system functions
addpath(genpath('SNS')); % include folder with inertial navigation system functions
addpath(genpath('Integrated')); % include folder with integrated navigation system functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));

filterTypeArray         = {'ukf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf'};

sampleTime              = 0.01; % seconds
simulationTime          = 1*0.05*60; % y hours * 60 min * 60 sec
simulationNumber        = round(simulationTime / sampleTime);
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;

accBiasMu               = zeros(3, 1);  % [km / sec^2]
accBiasSigma            = 5e-8*[1;1;1]; % [km / sec^2]
accNoiseVar             = 4e-8*[1;1;1]; % [km / sec^2]
accScale                = 1e-5*diag(3);

gyroBiasMu              = zeros(3, 1);  % [rad / sec]
gyroBiasSigma           = 5e-5*[1;1;1]; % [rad / sec]
gyroNoiseVar            = 4e-4*[1;1;1]; % [rad / sec]
gyroScale               = 1e-5*diag(3);

initialAcceleration     = zeros(3,1);
initialAngularVelocity  = zeros(3,1);
initialQuaternion       = [1;0;0;0];

%%  
accelerationInBodyFrame = AccelerationInBodyFrame(simulationNumber, ... simulationNumber
    sampleTime, ... sampleTime [sec]
    initialAcceleration, ... mu [km/sec^2]
    1.75*1e-7*[1; 1; 1] ... sigma [km/sec^2]
);

angularVelocityInBodyFrame = AngularVelocityInBodyFrame(simulationNumber, ... simulationNumber
    sampleTime, ... sampleTime
    initialAngularVelocity, ... mu
    1.2*1e-3*[1;1;1] ... sigma
);

tEpoch = 0.1465;
initialOrbit = loadInitialOrbit();

tic;

simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, tEpoch);
spaceshipStateTrue = simulator.Simulate(time, sampleTime, 1);
fprintf('Simulating true trajectory: ');

toc;


satellitePhaseState = SatellitePhaseSpace(initialOrbit, simulationNumber);
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
insInitArgs.T_till_current_epoch         = tEpoch;
insInitArgs.sampleTime                   = sampleTime;
insInitArgs.accNoiseVar                  = accNoiseVar;
insInitArgs.gyroNoiseVar                 = gyroNoiseVar;
insInitArgs.accScale                     = accScale;
insInitArgs.gyroScale                    = gyroScale;

snsInitArgs.Trajectory = spaceshipStateTrue.Trajectory';
snsInitArgs.Velocity   = spaceshipStateTrue.Velocity';
            
initCov = [(3)^2*eye(3) zeros(3, 19); ... distance error [km]^2
    zeros(3, 3) (5e-3)^2*eye(3) zeros(3, 16); ... velocity error [km/sec]^2
    zeros(4, 6) (1e-4)^2*eye(4) zeros(4, 12); ... quaternion error [-]
    zeros(3, 10) diag((insInitArgs.accBiasSigma.*insInitArgs.accBiasSigma)) zeros(3, 9); ... acceler bias [km/sec^2]^2
    zeros(3, 13) diag((insInitArgs.gyroBiasSigma.*insInitArgs.gyroBiasSigma)) zeros(3, 6); ... gyro bias [rad/sec]^2
    zeros(3, 16) 1e-10*eye(3) zeros(3, 3); ... acceler scale factor [-]
    zeros(3, 19) 1e-10*eye(3) ... gyro scale factor [-]
];
            
for j = 1:length(filterTypeArray)
    estimatorType = filterTypeArray(j);
    fprintf('estimator: %s\n', estimatorType{1});
    tic;
    
%     parfor k = 1:iterationNumber    
    for k = 1:iterationNumber    
        insSns = IntegratedInsSns(insInitArgs, snsInitArgs, initCov, initialAcceleration, initialAngularVelocity, initialQuaternion, time, sampleTime);        
        iterations(k) = insSns.Simulate(initialOrbit + [3*randn(3,1); 5e-3*randn(3,1); zeros(4,1)], estimatorType, 0, spaceshipStateTrue);
        
        fprintf ('iteration of %d: completed', k );
    end
       
    toc;    
end