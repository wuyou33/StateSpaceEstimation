close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('../../'));

set(0, 'defaultfigurecolor', [1 1 1]);
% ghqf not working due to dimension, too number of points :(
% {'ekf', 'ukf', 'cdkf', 'ckf', 'sckf', 'srukf', 'srcdkf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'pf', 'sppf', 'gspf', 'gmsppf'};

filterTypes  = {'srukf'};

date.day  = 17;
date.mon  = 11;
date.year = 2017;

initialOrbit = loadInitialOrbit();
m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);

tStart = '00:00:00.000';
tEnd = '00:20:00.000';

timeData = TimeExt(tStart, tEnd, 10, date, 1e12); % time data for integrated navigation system
timeDataSubSystem  = TimeExt(tStart, tEnd, 1, date, 1e12); % time data for inertial navigation system & satellite navigation system

iterationNumber             = 100;
drawIterations              = 0;
esitimatedParams            = 6;

mass                        = 200;
accBiasMu                   = zeros(3, 1);      % [km / sec^2]
accBiasSigma                = 5e-8*ones(3, 1);  % [km / sec^2]
accNoiseVar                 = 1e-6*ones(3, 1);  % [km / sec^2]
accScale                    = 5e-5*eye(3);      % [-]
gyroBiasMu                  = zeros(3, 1);      % [rad / sec]
gyroBiasSigma               = 5e-7*ones(3, 1);  % [rad / sec]
gyroNoiseVar                = 1e-6*ones(3, 1);  % [rad / sec]
gyroScale                   = 5e-5*eye(3);      % [-]
levelArm                    = zeros(3, 1);      % [g] where g is Earth gravity acc
angularAccelerBodyFrame     = zeros(3, 1);
gyroGSensitiveBias          = zeros(3);

initialAcceleration         = zeros(3, 1);          % [km/sec^2]
initialAngularVelocity      = zeros(3, 1);          % [rad/sec]
initialQuaternion           = initialOrbit(7:10);   % [-]

accelerationSigma           = 2e-5*ones(3, 1);      % [km/sec^2]
angularVelocitySigma        = 1e-7*ones(3, 1);      % [rad/sec]

insTrajInitErrorKm          = 5*ones(3, 1);         % [km]
insVelInitErrorKmSec        = 5e-3*ones(3, 1);      % [km/sec]
insQuaternionInitError      = 1e-9*ones(4, 1);      % [-]

accScaleFactorInitError     = 1e-5*ones(3, 1);      % [-]
gyroScaleFactorInitError    = 1e-5*ones(3, 1);      % [-]
gyroBiasSigmaInitError      = 1e-3*gyroBiasSigma;   % [rad / sec]
accBiasSigmaInitError       = 1e-3*accBiasSigma;    % [km / sec^2]
reconciliationTime          = 1e12;                 % [sec]

ansMeanTrajectory           = 0;                    % [km]
ansMeanVelocity             = 0;                    % [(km / sec)]

ansSigmaTrajectoryList = [(3.45)^2;    (1.1)^2;]; % [km^2]
ansSigmaVelocityList   = [(1.42e-3)^2; (1e-3)^2]; % [(km / sec)^2]
aTeffLegend    = {'AT = 5e3', 'AT = 1e5'};

insInitialState = initialOrbit + [insTrajInitErrorKm.*randn(3, 1); insVelInitErrorKmSec.*randn(3, 1); [1; 0; 0; 0] + insQuaternionInitError.*randn(4, 1)];
insInitialState(7:10) = quaternionNormalize(insInitialState(7:10));

errors = zeros(esitimatedParams, timeDataSubSystem.SimulationNumber);

%{
alphaList   = [(0.001)^2; (0.01)^2; (0.1)^2; (1)^2; (10)^2; (100)^2; (1000)^2];
betaList    = [(0.1)^2; (1)^2; (10)^2;];
gammaList   = [1];
%}

%%{
alphaList   = [1];
betaList    = [1];
gammaList   = [(0.001)^2; (0.01)^2; (0.1)^2; (1)^2; (10)^2; (100)^2; (1000)^2];
%}

initCov = [diag(insTrajInitErrorKm).^2 zeros(3, 19); ... distance error [km]^2
    zeros(3, 3) diag(insVelInitErrorKmSec).^2 zeros(3, 16); ... velocity error [km/sec]^2
    zeros(4, 6) diag(insQuaternionInitError).^2 zeros(4, 12); ... quaternion error [-] 2.--- 5e-3
    zeros(3, 10) diag(accBiasSigmaInitError).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2 --- 5e-2
    zeros(3, 13) diag(gyroBiasSigmaInitError).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
    zeros(3, 16) diag(accScaleFactorInitError).^2 zeros(3, 3); ... acceler scale factor [-]
    zeros(3, 19) diag(gyroScaleFactorInitError).^2 ... gyro scale factor [-]
    ];

tic;

estimatorType = filterTypes(1);
iinsProcCov = getIntegratedInsProcCov(estimatorType, accBiasSigma, gyroBiasSigma);

iinsInitArgs.initialParams(1:3)           = accBiasMu;
iinsInitArgs.initialParams(4:6)           = accBiasSigma;
iinsInitArgs.initialParams(7:9)           = gyroBiasMu;
iinsInitArgs.initialParams(10:12)         = gyroBiasSigma;
iinsInitArgs.initialParams(13:15)         = initialAcceleration;
iinsInitArgs.initialParams(16:18)         = initialAngularVelocity;
iinsInitArgs.initialParams(19:28)         = insInitialState;
iinsInitArgs.initialParams(29)            = timeData.SampleTime;
iinsInitArgs.initialParams(30)            = timeData.StartSecond;
iinsInitArgs.processNoiseMean             = zeros(22, 1);
iinsInitArgs.processNoiseMean(11:16)      = [accBiasMu; gyroBiasMu];
iinsInitArgs.observationNoiseMean         = [zeros(3, 1); zeros(3, 1)];

iterations = zeros(iterationNumber, esitimatedParams, timeDataSubSystem.SimulationNumber);

for at = 1:length(aTeffLegend)
    for b = 1:length(betaList)
        for a = 1:length(alphaList)
            for g = 1:length(gammaList)
                ansSigmaTrajectory = ansSigmaTrajectoryList(at);
                ansSigmaVelocity   = ansSigmaVelocityList(at);
                alpha   = alphaList(a);
                beta    = betaList(b);
                gamma   = gammaList(g);
                aTeffLabel = aTeffLegend(at);
                
                iinsInitArgs.processNoiseCovariance       = beta*iinsProcCov;
                iinsInitArgs.observationNoiseCovariance   = gamma*[[ansSigmaTrajectory*eye(3) zeros(3)]; [zeros(3) ansSigmaVelocity*eye(3)]];
                
                gravModel = m_fitSolarSystemGravityModel(timeDataSubSystem.SampleTime, timeDataSubSystem.SimulationNumber);
                parfor j = 1:iterationNumber
                    accelerationInBodyFrame     = AccelerationInBodyFrame(timeDataSubSystem, initialAcceleration, accelerationSigma);
                    angularVelocityInBodyFrame  = AngularVelocityInBodyFrame(timeDataSubSystem, initialAngularVelocity, angularVelocitySigma);
                    
                    simulator = TrajectoryPhaseSpaceSatelliteSimulator(initialOrbit, accelerationInBodyFrame, angularVelocityInBodyFrame, timeDataSubSystem, mass);
                    xTrueState = simulator.simulate(iterationNumber == 1);
                    insInitArgs = struct('accBiasMu', accBiasMu, 'accBiasSigma', accBiasSigma, 'gyroBiasMu', gyroBiasMu, 'gyroBiasSigma', gyroBiasSigma, 'visualize', 0,...
                        'timeData', timeDataSubSystem, 'accNoiseVar', accNoiseVar, 'gyroNoiseVar', gyroNoiseVar, 'accScale', accScale, 'gyroScale', gyroScale,...
                        'levelArm', levelArm, 'angularAccelerBodyFrame', angularAccelerBodyFrame, 'gyroGSensitiveBias', gyroGSensitiveBias,...
                        'gravityModel', gravModel,'mass', mass, 'accelerationInBodyFrame', accelerationInBodyFrame, 'angularVelocityInBodyFrame', angularVelocityInBodyFrame);
                    
                    trueState = xTrueState;
                    insSnsInitState = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]; [0; 0; 0]; 5e-6*ones(3, 1); 5e-6*ones(3, 1)];
                    
                    % gaussian | markov1 | wiener | exp | realImitator
                    fns = NavImitator(trueState.State, ansSigmaTrajectory, ansMeanTrajectory, ansSigmaVelocity, ansMeanVelocity, 'gaussian');
                    ins = initInertialNavigationSystem('init', insInitArgs);
                    
                    iins = IntegratedIns(ins, fns, timeData, iinsInitArgs, timeDataSubSystem, reconciliationTime);
                    iinsState = iins.evaluate(insSnsInitState, alpha*initCov, insInitialState, estimatorType, iterationNumber == 1, iterationNumber == 1);
                    
                    quatErr = trueState.Rotation - iinsState.Rotation;
                    errTraj = vectNormError(trueState.Trajectory, iinsState.Trajectory, 1e3);
                    errVel  = vectNormError(trueState.Velocity, iinsState.Velocity, 1e3);
                    
                    iterations(j, :, :) = [errTraj; errVel; quatErr];
                end
                
                if iterationNumber > 1
                    for i = 1:timeDataSubSystem.SimulationNumber
                        errors(:, i) = ( sum( ( squeeze(iterations(:, :, i))' ).^2, 2) / iterationNumber ).^0.5;
                    end
                end
                
                fprintf('%s; a = %d; b = %d; g = %d\n', aTeffLabel{1}, sqrt(alpha), sqrt(beta), sqrt(gamma));
                fprintf('\tRMS trajectory: %d\n', errors(1, end));
                fprintf('\tRMS velocity: %d\n', errors(2, end));
            end
        end
    end
end
toc;
fprintf('DONE\n\n');
