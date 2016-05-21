close all force; clc;
clear all;
addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('TOAMeasurementUnit')); % include folder with Pulsar & Quasar navigation
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));
addpath(genpath('Ephemeris'));
addpath(genpath('XNAV'));

sampleTime              = 0.0001; % seconds
simulationNumber        = fix(1*0.5/sampleTime); % x minutes * y seconds / sampleTime
simulationTime          = sampleTime*simulationNumber;
time                    = (1:simulationNumber) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;

filterTypeArray         = {'srcdkf'}; %{'ukf', 'cdkf', 'ckf', 'sckf', 'srukf','srcdkf', 'pf'};

T_till_current_epoch = 0.1465;
xRaySourceCount      = 7;
backgroundPhotnRate  = 5.9;
timeBucket           = 1e5; % 1e5
detectorArea         = 1;

inter = loadEphemeris('earth', simulationNumber, 60/sampleTime);
% orig = loadEphemeris('earth', simulationNumber);
t = ((inter.x(simulationNumber) - inter.x(simulationNumber-1)) / (sampleTime)) - inter.vx(simulationNumber);
display(t)
