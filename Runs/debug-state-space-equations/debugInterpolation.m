close all; clc; clearvars; clear memoize; % clear memoize required for memoization

addpath(genpath('./'));

matrix = randn(2);
qr(matrix, 0)
chol(matrix)

% return
sampleTime              = 0.0001; % seconds
simulationNumber        = fix(1*0.5/sampleTime); % x minutes * y seconds / sampleTime
simulationTime          = sampleTime*simulationNumber;
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;

inter = loadEphemeris('earth', simulationNumber, 60/sampleTime);
% orig = loadEphemeris('earth', simulationNumber);
interpError = ((inter.x(2) - inter.x(1)) / (sampleTime)) - inter.vx(2);
display(interpError)

interpErrorEnd = ((inter.x(simulationNumber) - inter.x(simulationNumber-1)) / (sampleTime)) - inter.vx(simulationNumber);

display(interpErrorEnd);
