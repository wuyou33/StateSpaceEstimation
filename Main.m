% clear all; 
close all force; clc;

% 
% myfunction = @(x, y) [...
%     x^2;
%     x*y...
%     ];

addpath(genpath('func_Vova') ); % ���������� ����� � ��������� �� ����
addpath(genpath('My_func_aux') ); % ���������� ����� � ����� ���������
addpath(genpath('My_func_general') ); % ���������� ����� � ����� ���������
addpath(genpath('OrbitalMotion')); % include folder with orbital motion functions
addpath(genpath('InertialMeasurementUnit')); % include folder with inertial measurement functions
addpath(genpath('InertialNavigationSystem')); % include folder with inertial navigation system functions
addpath(genpath('SNS')); % include folder with inertial navigation system functions
addpath(genpath('Integrated')); % include folder with integrated navigation system functions
addpath(genpath('Statistics')); % include folder with Statistics functions
addpath(genpath('TOAMeasurementUnit')); % include folder with Pulsar & Quasar navigation
addpath(genpath('StateSpaceEstimation')); % include folder with State Space Estimation algorithms
addpath(genpath('Utils'));

filterTypeArray         = {'ukf'}; %{'ukf', 'ckf', 'srckf', 'cdkf','srukf','srcdkf'};
sampleTime              = 1; % seconds
simulationTime          = 1*60*60; % y hours * 60 min * 60 sec
simulationNumber        = round(simulationTime / sampleTime);
time                    = (1: simulationNumber +1) * sampleTime;
timeMinutes             = time / 60;
iterationNumber         = 1;
accBiasMu               = 1e-14*[1;1;1];
accBiasSigma            = 5e-11*[1;1;1];
gyroBiasMu              = 1e-14*[1;1;1];
gyroBiasSigma           = 5e-11*[1;1;1]; 
initialAcceleration     = zeros(3,1);
initialAngularVelocity  = zeros(3,1);
initialQuaternion       = [1;0;0;0];

%%  
accelerationInBodyFrame = AccelerationInBodyFrame([0; 0; 0], ... mu
    1.2*1e-7*randn(3,1), ... sigma    
    simulationNumber ... simulationNumber
);

angularVelocityInBodyFrame = AngularVelocityInBodyFrame(1e-5*randn(3, 1), ... mu
    1.2*1e-3*randn(3,1), ... sigma    
    simulationNumber ... simulationNumber
);

T_till_current_epoch = 0.1465;
% initial = [-21223.9926714100; -42868.3395565589; 0; 2.58697266398439; % -1.28080278894642; 0; 1; 0; 0; 0]; % GEO sat
initial = [-6158.34755458333; -4063.45976435815; 0; 3.98653859590107; -6.04177022463943; 1.27633791914791; 1; 0; 0; 0]; % LO sat

tic;
simulator = TrajectoryPhaseSpaceSatelliteSimulator(initial, ...
    accelerationInBodyFrame, ... accelerationInBodyFrame
    angularVelocityInBodyFrame, ... angularVelocityInBodyFrame
    T_till_current_epoch);

satellitePhaseStateTrue = simulator.Simulate(time, sampleTime);

% satellitePhaseStateTrue = EquationOfMotionSolver(y0_sc_true_initial, accelerationInBodyFrame.Acceleration, angularVelocityInBodyFrame.Velocity, time, T_till_current_epoch, sampleTime);
fprintf('Simulating true trajectory: ');
toc;

SatelliteOrbitVisualization(satellitePhaseStateTrue);

%% simulate INS
satellitePhaseState = SatellitePhaseSpace(initial, simulationNumber + 1);
iterations = repmat(satellitePhaseState, iterationNumber, 1);

insInitArgs.accBiasMu                    = accBiasMu;
insInitArgs.accBiasSigma                 = accBiasSigma;
insInitArgs.gyroBiasMu                   = gyroBiasMu;
insInitArgs.gyroBiasSigma                = gyroBiasSigma;
insInitArgs.accelerationInBodyFrame      = accelerationInBodyFrame;
insInitArgs.angularVelocityInBodyFrame   = angularVelocityInBodyFrame;
insInitArgs.simulationNumber             = simulationNumber;
insInitArgs.timeMinutes                  = timeMinutes;
insInitArgs.visualize                    = 1;
insInitArgs.T_till_current_epoch         = T_till_current_epoch;

snsInitArgs.Trajectory = satellitePhaseStateTrue.Trajectory';
snsInitArgs.Velocity   = satellitePhaseStateTrue.Velocity';

initialCov = [1e2*eye(3) zeros(3, 19); ... distance error
    zeros(3, 3) 1e-2*eye(3) zeros(3, 16); ... velocity error
    zeros(4, 6) 1e-4*eye(4) zeros(4, 12); ... quaterni error
    zeros(3, 10) diag((accBiasSigma.*accBiasSigma)) zeros(3, 9); ... acceler bias
    zeros(3, 13) diag((gyroBiasSigma.*gyroBiasSigma)) zeros(3, 6); ... gyro bias
    zeros(3, 16) 1e-6*eye(3) zeros(3, 3); ... acceler scale factor
    zeros(3, 19) 1e-6*eye(3) ... acceler scale factor
];

initialState = sqrt(initialCov)*randn(22, 1);

tic;

% parfor k = 1:iterationNumber
for j = 1:length(filterTypeArray)
    estimatorType = filterTypeArray(j);
    fprintf('estimator: %s\n', estimatorType{1});
    
    for k = 1:iterationNumber    
        copyTime = time;    
        ins = initInertialNavigationSystem('init', insInitArgs);
        sns = initSatelliteNavigationSystem('init', snsInitArgs);

        gssmInsSnsInitArgs.initialParams                 = [accBiasMu; accBiasSigma; gyroBiasMu; gyroBiasSigma; initialAcceleration; initialAngularVelocity; initialQuaternion, sampleTime, time(1)];
        gssmInsSnsInitArgs.processNoiseMean              = [accBiasMu; gyroBiasMu];
        gssmInsSnsInitArgs.processNoiseCovariance        = [diag(accBiasSigma.*accBiasSigma) zeros(3,3); zeros(3,3) diag(gyroBiasSigma.*gyroBiasSigma)];
        gssmInsSnsInitArgs.observationNoiseMean          = zeros(6, 1); 
        gssmInsSnsInitArgs.observationNoiseCovariance    = [10*10*eye(3) zeros(3,3); zeros(3,3) .01*.01*eye(3)];

        insErrorDynamicStateSpaceModel = gssmInsSns('init', gssmInsSnsInitArgs);
        args.type  = 'state';
        args.tag   = 'State estimation for loosely coupled Ins & Sns integrated system';
        args.model = insErrorDynamicStateSpaceModel;

        % todo: debug inference noise generator
        [processNoise, observationNoise, inferenceDataSet] = inferenceNoiseGenerator(inferenceDataGenerator(args), estimatorType);
        
        switch estimatorType{1}
            case 'ukf'
                alpha = 1; % scale factor (UKF parameter)
                beta  = 2; % optimal setting for Gaussian priors (UKF parameter)
                kappa = 0; % optimal for state dimension=2 (UKF parameter)
                
                inferenceDataSet.spkfParams = [alpha beta kappa];
                
            case 'srukf'
                alpha = 1; % scale factor (UKF parameter)
                beta  = 2; % optimal setting for Gaussian priors (UKF parameter)
                kappa = 0; % optimal for state dimension=2 (UKF parameter)
                
                inferenceDataSet.spkfParams = [alpha beta kappa];
                Sx = chol(Px)';
                
            case 'cdkf'
                inferenceDataSet.spkfParams = sqrt(3); % scale factor (CDKF parameter h)
                
            case 'srcdkf'
                alpha = 1;         % scale factor (UKF parameter)
                beta  = 2;         % optimal setting for Gaussian priors (UKF parameter)
                kappa = 0;         % optimal for state dimension=2 (UKF parameter)
                
                Sx = chol(Px)';
                
                inferenceDataSet.spkfParams = [alpha beta kappa];
            otherwise
                % do nothing by default
        end
    %     tic;
        state = initialState;
        covState = initialCov;
        insMeasurement = initial;
        
        for i = 2:1:simulationNumber + 1          
            insMeasurement = ins.Simulate(insMeasurement, i, sampleTime, copyTime(i));
            snsMeasurement = sns.Simulate(i);
            observation = insMeasurement(1:6) - snsMeasurement;
            
            modelParams.params(1:3)    = inferenceDataSet.model.params(1:3);   % accelerationBiasMu
            modelParams.params(4:6)    = inferenceDataSet.model.params(4:6);   % accelerationBiasSigma
            modelParams.params(7:9)    = inferenceDataSet.model.params(7:9);   % gyroBiasMu
            modelParams.params(10:12)  = inferenceDataSet.model.params(10:12); % gyroBiasSigma
            modelParams.params(13:15)  = ins.GetAcceleration(i); 
            modelParams.params(16:18)  = ins.GetAngularValocity(i);
            modelParams.params(19:22)  = insMeasurement(7:10);                 % quaternion
            modelParams.params(23)     = inferenceDataSet.model.params(23);    % sampleTime
            modelParams.params(24)     = time(i);
            
            inferenceDataSet.model.setParams(modelParams);
            
            [state, covState, processNoise, observationNoise, ~] = ukf(state, covState, processNoise, observationNoise, observation, inferenceDataSet);
            iterations(k).AddPhaseState(state, i);
            % TODO: check it
            insMeasurement = insCorrection(insMeasurement, state(1:10));
        end

        SatelliteOrbitVisualization(iterations(k));
        fprintf ('thread of %d: ', k );    
    %     toc;
    end

    fprintf('Simulation: ');
    toc;

    tic;
    errorReducer = ErrorReducer(iterations, satellitePhaseStateTrue, 10, simulationNumber);
    rmsd = errorReducer.RMSD();
    display('RMSD calculation: ');
    toc;

    figure(); 
        plot(timeMinutes(2:end)', [sqrt(rmsd(:, 1))'; sqrt(rmsd(:, 2))'; sqrt(rmsd(:, 3))']); 
        title('Root-mean-square deviation displacement'); 
        ylabel('RMSD displacement, m');
        xlabel('time, min');
        legend('x axis', 'y axis', 'z axis');
        grid on;

    figure(); 
        plot(timeMinutes(2:end)', [sqrt(rmsd(:, 4))'; sqrt(rmsd(:, 5))'; sqrt(rmsd(:, 6))']); 
        title('Root-mean-square deviation velocity'); 
        ylabel('RMSD velocity, m/s');
        xlabel('time, min');
        legend('x axis', 'y axis', 'z axis');
        grid on;   
end
%%

%{
%% =========================== ������������� ���������� ������======================
%��������� ������ ��� ������������ ���� � �� �� ���������� �����, ����� ���������� ��������
s=RandStream('mt19937ar','Seed',785);%784
RandStream.setGlobalStream(s);
% RandStream.getGlobalStream


%% =========================== ������������� ��������� ���======================
%���� ����� ��������� ������� ���
fileAlmanah=[cd,'\data\MCCT_140825.AGL'];   %���� � ����� ���������
%������������� ���������
alm_gln = almanah(fileAlmanah); % �������� ���������, ���������� �������� �������� ��� �����-�� ����

alm_testSetilate = almanahReciver(1);%�������� ������ � �������

%----------------------------------------------------------------
%% ========================= ����� ��������� ==================================
global GL_Hmin_ion%[��] ������ ���� ������� ���������,��� �� �� ����� 
GL_Hmin_ion=500;%[��] ������ ���� ������� ���������,��� �� �� ����� 

global GL_Ae %[��] - �������������� ������ �����
GL_Ae=6378.136;%[��] - �������������� ������ �����

Year_ref=1996;%������� ���, �� �������� ������������� N4 (�� ���)
global GL_W_rot_Earth %[���/�] - ������� �������� �������� �����
GL_W_rot_Earth=7.2921151467e-5;%[���/�] - ������� �������� �������� �����

global GL_MuE%[��^3/�^2] - ��������� ��������������� ���� �����
GL_MuE=398600.4418;%[��^3/�^2] - ��������� ��������������� ���� �����

global J2_0;% - ������ ��������� ��������� �������������
J2_0=1082625.75e-9;% - ������ ��������� ��������� �������������

global GL_Temp;%[K] - ����������� ��������
GL_Temp=180;%[K] - ����������� ��������

global GL_T_Navig_message
GL_T_Navig_message=1e-3;%[��] - ����� �������������� 

global GL_T_chip% [���]  - ������������ 1 ���� ��� ������������� ����
N_chip=511;%[��] - ����� ����� ���
GL_T_chip=GL_T_Navig_message/N_chip;% [���]  - ������������ 1 ���� ��� ������������� ����

global GL_Fdop_cell%[��] - ������ ���������� �� ������� ������� � ���������� ������ � �������
GL_Fdop_cell=500;%[��] - ������ ���������� �� ������� ������� � ���������� ������ � �������

global GL_DN_GLN %[��] ��������� �������������� ��� �������

GL_DN_GLN = [0.0    11.0; 4.0    11.5; 8.0    13.0; 10.0   13.6; 12.0   14.0; 16.0   13.2; 20.0   10.8;... %new
    24.0    6.0; 28.0   -1.0; 30.0   -3.0; 36.0   -6.0; 44.0   -8.0; 50.0  -10.0; 52.0  -12.0; 54.0  -15.0; 60.0  -11.0; ...
    64.0   -7.0; 70.0   -6.0; 80.0   -7.0; 90.0   -10.0];

global GL_DN_ARN %[dB] ��������� ������������� ��

GL_DN_ARN = [0	14.8; 5	14.5; 10	14.0; 15	13; 20	10.5; 25	8; 30	3; 35	-6; 40	-6; 45	-2; 50	0; 55	-3;...%new
    60	-6; 65	-18; 70	-12;75	-10; 80	-12; 85	-18; 90	-30];

global GL_DN_ARN_LO %[dB] ��������� ������������� �� ��� ��
GL_DN_ARN_LO = [90 3;100 3;110 3;120 3;130 3;140 3;150 3; 160 3;170 3;180 3];% 190	3; 200	3; 210	3; 220	3;230 3;240	3;250 3;260	3; 270	3];%[dB] ��������� ������������� �� �� �� ��� � �������


global GL_Threshold;%[dB-Hz] - ����� ��� ������� ����� ���������, ��� ������� ��������� ��� ���
GL_Threshold=33;%[dB-Hz] - ����� ��� ������� ����� ���������, ��� ������� ��������� ��� ��� !!!!!!!!!!!!

%================================
%% C_N0=GL_Threshold%40;%[��-��] -             ���������� �� ������� !!!!!!!!!!!!!!!!!!!!
% C_N0=40;%[��-��] -                                      ���������� �� ������� !!!!!!!!!!!!!!!!!!!!
C_N0=zeros(24,1);%[��-��] -                                      ���������� �� ������� !!!!!!!!!!!!!!!!!!!!
C=40;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_stat=1;%����� ���������� ��� ������ ���������� !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

time_refresh_data=15*60;%[���] = 15��� - �����, ����� ������� ���������� ��������� ������ � ������������ ���� � ������
time_sat_change=60;%[���] �����, ����� ������� �� ������ ���

date_formatIn = 'dd.mm.yyyy'; %������ ����� ������ � ����
time_formatIn = 'HH:MM:SS.FFF'; %������ ����� ������ � �������

 %!!!!!!!!!!!!!!!!!!!!!!!!!!!
% T_dT_corr=1;%[���] - �������� �������, ������ �������� ����� ����������� ���������� ��������� ��=sc
T_dT_sc=1/400*1e0;%[���] - �������� ������������� �� �������, � ������� �� ����� ������ �� ��� ���������� ��������� ���=nsc. ��� �� ��� �������� �������, � ������� ���� ����������!!!
T_accum=1e-3;%T_dT_sc;%����� �������������
T_dT_U=15000*1e-3;%[���] - ����� ��������� ����������

Vmax_err=0.04;%[��/�] - ��� ������ ��������� ���� ���������� �� ��������
Rmax_err=0.25;%[��] - ��� ������ ��������� ���� ���������� �� ���������

%% ===================== ����� ��������� ������� ================================

time_ref='14:30:00.000';%����� ������ ����������
time_end_eval='14:35:00.00';%����� � ����ר���� ����� (����������!!!!! ���)

date_start_eval_mass=datevec(time_ref,time_formatIn);%������ � ��������, �� ������� ��� ���� ��������� 

fields_date=fieldnames(alm_gln.date);%��������,����� ���� ���� � �������� � ����� �� ���������
for i=1:size(fields_date,1)
    date_start_eval_mass(4-i)=alm_gln.date.( char( fields_date(i) ) ); % ������������ ��������� ������ � ����� � ��������� �������� �������
end;

date_end_eval =[date_start_eval_mass(1:2) date_start_eval_mass(3)+0];%���� � ����ר���� ����� (����������!!!!! ���) %!!! ������ ��� �������������. �������� 1 � 3�� �����-����������� �� 1 �����

Date_end_eval_mass=datevec(time_end_eval,time_formatIn);%������ � �������� � ����ר���� �����
Date_end_eval_mass(1:3)=date_end_eval(:);%�������� ������ � ����� � �������� � ����ר���� �����

clear i
%----------------------------------------------------------------
%% ====================== ������ ��������� ����������===========================

N4=fix( (date_start_eval_mass(1)-Year_ref)/4 )+1;%N4 - ����� ��������� �� 4 ����(����� ������� �����������) � ������� ��������� ������� ����
%� ������� datenum ������ ������: yyyy.mm.dd
Nt=datenum( date_start_eval_mass(1:3) ) - datenum([Year_ref+(N4-1)*4 1 1])+1;%Nt - ����� ����� ������ 4�������� ����������� �������

% [Day_cur, Month_cur, Year_cur]=N4Nt2currentDate(N4,Nt);%�������� ���� �� ������� N4 � Nt � ������ dd.mm.yyyy

JD=JD_count_fun(N4,Nt);%��������� ������� ��������� ����

Delta_dateNtime=Date_end_eval_mass-date_start_eval_mass;%����������� ������� �� ������� �/� ����� � �������� ����� ���������� � ��������� ����� � ��������

T_Ttotal_eval=Delta_dateNtime(3)*3600*24+T_fun( Delta_dateNtime(4:6) );%[���] - ����� ������ � ���� �������

T_Tstart=0*T_fun( date_start_eval_mass(4:6) );%[���] - ��������� ����� �������� ��������(������ ������ ����������)
T_Tend=T_Ttotal_eval+T_Tstart;


clear Delta_dateNtime

%----------------------------------------------------------------
%% ==================== ���������� ��������� ��� � �� �� ������� ���������=====================

nsmax = 24;           %���-�� ������������ ���������

%=========������ ������� ���������� �� � ��� � �������� �� �� ����������� �� � ������������ ��============

tic
[Satpos_xyz_Rec_mean,Satpos_xyz_gln,Sun,Moon] = SC_and_NSC_trajectory_calc(...
    T_Ttotal_eval,...
    T_dT_sc, ...
    time_refresh_data, ...
    JD, ...
    T_Tstart, ...
    T_Tend, ...
    alm_testSetilate, ...
    date_start_eval_mass, ...
    Nt,alm_gln);
toc
figure(); plot(Satpos_xyz_Rec_mean.x);
figure(); plot(Satpos_xyz_Rec_mean.y);
figure(); plot(Satpos_xyz_Rec_mean.z);
figure(); plot(Satpos_xyz_Rec_mean.vx);
figure(); plot(Satpos_xyz_Rec_mean.vy);
figure(); plot(Satpos_xyz_Rec_mean.vz);
figure(); plot(Satpos_xyz_Rec_mean.q);
figure(); plot(Satpos_xyz_Rec_mean.qi);
figure(); plot(Satpos_xyz_Rec_mean.qj);
figure(); plot(Satpos_xyz_Rec_mean.qk);
return
%{
%----------------------------------------------------------------------
%% =======����������� ���������� �� � ������� ���������� ������������ �� �������======================
%������ ��� ������ ����������
N_max=T_Ttotal_eval/T_dT_sc+1;% ������ ����� ����� ���������� �����
Satpos_xyz_Rec_current(1:N_stat)=struct('x',zeros(1,N_max),'y',zeros(1,N_max),'z',zeros(1,N_max),'vx',zeros(1,N_max),'vy',zeros(1,N_max),'vz',zeros(1,N_max));
State_space_error=zeros(6,N_stat);

% matlabpool open %��� ������������ ����������

tic
for i=1:N_stat
    State_space_error(:,i)=1*[Rmax_err*randn; Rmax_err*randn; Rmax_err*randn; Vmax_err*randn; Vmax_err*randn; Vmax_err*randn];%������ �� ����������� � ��������� � ��������� ������ �������!!!!!!!!!!!!!!!
%     State_space_error(:,i)=[-0.379039417699431;0.213644210733933;0.390210249361874;-0.0220798398887709;0.0939310296257019;0.00691360219067245];
    
    while sqrt( sum( ( ( State_space_error(4:6,i) ).*1000).^2 ) )>120;
        State_space_error(4:6,i)=[Vmax_err*randn; Vmax_err*randn; Vmax_err*randn];%������ �� ����������� � ��������� � ��������� ������ �������!!!!!!!!!!!!!!!
    end;
    
    y0_sc_current=[Satpos_xyz_Rec_mean.x(1);Satpos_xyz_Rec_mean.y(1);Satpos_xyz_Rec_mean.z(1);...
                             Satpos_xyz_Rec_mean.vx(1);Satpos_xyz_Rec_mean.vy(1);Satpos_xyz_Rec_mean.vz(1)]+...
                             State_space_error(:,i);% ����� ��������� ������� ��� ��  

   Satpos_xyz_Rec_current(i)=SC_current_trajectory_calc(T_Ttotal_eval,time_refresh_data,T_dT_sc,T_Tstart,T_Tend,...
        Sun,Moon,y0_sc_current);                      
    
end;
toc
% matlabpool close %��� ������������ ����������

% r_sc_current=sqrt( Satpos_xyz_Rec_current.x(:).^2+Satpos_xyz_Rec_current.y(:).^2+Satpos_xyz_Rec_current.z(:).^2 )-GL_Ae;
% figure;
% plot(r_sc_current);

% r_sc_mean=sqrt( Satpos_xyz_Rec_mean.x(:).^2+Satpos_xyz_Rec_mean.y(:).^2+Satpos_xyz_Rec_mean.z(:).^2 )-GL_Ae;
% figure;
% plot(r_sc_mean);


%{
% R=zeros(1,N_max);

Err_coord=Coord_Inert_sc_true-Coord_Inert_sc_current;%��������� ������� ������� ���������� �� ��������
Err_vel=Vel_Inert_sc_true-Vel_Inert_sc_current;%��������� ������� ������� �������� �� ��������

r_sc_true=sqrt( Coord_Inert_sc_true(1,:).^2+Coord_Inert_sc_true(2,:).^2+Coord_Inert_sc_true(3,:).^2 )-GL_Ae;
r_sc_current=sqrt( Coord_Inert_sc_current(1,:).^2+Coord_Inert_sc_current(2,:).^2+Coord_Inert_sc_current(3,:).^2 )-GL_Ae;
r_nsc=sqrt( Coord_Inert_nsc_gln(1,:).^2+Coord_Inert_nsc_gln(2,:).^2+Coord_Inert_nsc_gln(3,:).^2 )-GL_Ae;
%}

clear Vel_Inert_sc_current Coord_Inert_sc_current Vel_Inert_sc_true Coord_Inert_sc_true ti_current i K k Y 
clear Sun Moon start stop Vel_Inert_nsc_gln Coord_Inert_nsc_gln y0_sc ns S date_current_time_mass r_nsc
clear hour minute sec month ERA GMST j y0_sc_true time_formatIn flag r_sc_true r_sc_current
clear Err_vel Err_coord Ti_current
clear y0_sc_current
%% ========���������� ��������� ���������, �������� ������� � ������������� ������================


   
Out_param_current(1:nsmax,N_stat)=struct('vis',[],'Power',[],'Tz',[],'N_navig',[],'N_chip_Tau',[],'T_chip_resid',[],'T_navig_resid',[],'Fdop_1',[],'N_chip_Fdop',[],'Phi_HF_1',[],'cH',[],'Litera',[],'delta_fi0',[],'d_fd_dt',[]);
Out_param_mean(1:nsmax,1)=struct('vis',[],'Power',[],'Tz',[],'N_navig',[],'N_chip_Tau',[],'T_chip_resid',[],'T_navig_resid',[],'Fdop_1',[],'N_chip_Fdop',[],'Phi_HF_1',[],'cH',[],'Litera',[],'delta_fi0',[],'d_fd_dt',[]);
Visible_satellites=zeros(nsmax,N_max,N_stat);
Visible_satellites_reserv=zeros(nsmax,N_max,N_stat);


%����������� �������� ��������� ������� - ��������, ������, ���������, �����
%����� ����� � ���� ��������

tic
for i=1:N_stat
    v_fd=zeros(nsmax,N_max);
    for ns = 1:nsmax
        
        [Out_param_current(ns,i)] = Calc_param(Satpos_xyz_Rec_current(i), Satpos_xyz_gln(ns), alm_gln, ns);%�-��� ���������� ����������
%         [Out_param_current(ns,i)] = Calc_param(Satpos_xyz_Rec_mean, Satpos_xyz_gln(ns), alm_gln, ns);%�-��� ���������� ����������
              
        Visible_satellites(ns,:,i)=Out_param_current(ns,i).vis;%������� ������ ������� ���������
        
        delta_fi0=zeros(1,N_max);
        v_fd(ns,1:N_max-1)=diff(Out_param_current(ns,i).Fdop_1)/T_dT_sc;%[��/�]
        v_fd(ns,N_max)=v_fd(ns,N_max-1);
        Fdop_1=Out_param_current(ns,i).Fdop_1(:);
        
        delta_fi0(1)=2*pi*rand;%pi/2+10/180*pi*randn(1);%+T_dT_sc.*Fdop_1(1)*2*pi;
        for jj=2:N_max%��������� ���� � �����������
          delta_fi0(jj)=delta_fi0(jj-1)+T_dT_sc.*Fdop_1(jj-1)*2*pi+(T_dT_sc.^2./2).*v_fd(ns,jj-1)*2*pi; 
        end;      
        
        
        Out_param_current(ns,i).d_fd_dt=v_fd(ns,:);
        Out_param_current(ns,i).delta_fi0=delta_fi0;
          
%         rrrr(1)=Out_param_current(ns,i).delta_fi0(1);
%         rrrr(2:N_max)=Out_param_current(ns,i).delta_fi0(1)+time_mass.*Out_param_current(ns,i).Fdop_1(1:end-1)*2*pi+...
%             (time_mass.^2./2).*v_fd*2*pi;
          
    end;
end
toc


for ns = 1:nsmax
    [Out_param_mean(ns,1)] = Calc_param(Satpos_xyz_Rec_mean, Satpos_xyz_gln(ns), alm_gln, ns);%�-��� ���������� ����������
end;

%-----------------------------------
%������� ��� � ����������� ��������

tic
for i=1:N_stat
    
    [Visible_satellites(:,:,i),Out_param_current(:,i)]=double_NSC_discard_alg(nsmax,alm_gln,Visible_satellites(:,:,i),Out_param_current(:,i));    
    
end;
toc

%--------------------------------------------------
%������ 1 ��� �������� �� 3 ��� �� ������� ��� �������� ����������������� �� 3 ���
%{
visual=0;%������ ��� ��� ������ ��� �������� GDOP

G=zeros(N_stat,N_max);
flag=1;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
for i=1:N_stat
    while flag
        %������� 4 ��������� ���
        Visible_satellites_reserv(:,:,i)=NKA_selectioon(Visible_satellites(:,:,i),N_max,T_dT_sc,time_sat_change);
        for ns=1:nsmax
            Out_param_current(ns,i).vis(1:end)=Visible_satellites_reserv(ns,:,i);    
        end;
        G(i,:)=GDOP_new(Out_param_current(:,i),T_dT_sc,visual);
        if isempty( find(G(i,:)>20,1) )%���������, ������ �� �������������� ������, ���� ��, �� ������� ������
            flag=0;
            Visible_satellites(:,:,i)=Visible_satellites_reserv(:,:,i);
        end;
    end;
    flag=1;
end;
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');
%}
%--------------------------------------------------


%
visual=1;%������ ��� ��� ������ ��� �������� GDOP

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
tic
GDOP_new(Out_param_current(:,1),T_dT_sc,visual);
toc
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');


clear visual Err_Tau0 Err_Fdop0 K_new ns N_int_sec observ_quant num_sat num_sat2 select_sat_num combinations num_sat n_vis Visible_satellites_reserv 
clear v_fd
%}
%% ============================���������� �����=================================
%���������� �����
%{
T1=1;%��������� ������ ��� ����������
T2=10*1000;%�������� ������ ��� ����������
t_refresh=10;%����� ��������, ����� ������� ����� ��������� ����� ������� ���
satellite_orbit_visualization(T1,T2,t_refresh,Satpos_xyz_Rec_mean,Satpos_xyz_gln,Visible_satellites,nsmax)%��������� ���������� �����
%}



%{
����������:

1. ��� ������� v5.1 2008 - ���������� (2013-14��)
2. �.�. �������, "�������������� ��������� ���� ��� �� ����������", �:����� � �����, 1983 - �� ������������
3. Rudolph van der Merwe, "Sigma-Point Kalman Filters for Probabilistic Inference in Dynamic State-Space Models" PhD Thesis, 2004
4. �.�. ��������, "�������������� ������ ���������������� ������", 2007 � 
5. Chun Yang, "GPS Signal tracking with Kalman Filter based on joint code delay and carrier phase and frequency error discriminator" , ION 60th meeting 2004
6. �.�. �����, �.�. �������, "������� �������� ���������� � ����������������", ��� 4��, �:"������������", 2010.�
7. �.�. ���������, "�������� � ������ ������ ������������� ��������� �����", �, 1965  �
8. ������������� �����: http://en.wikipedia.org/wiki/Rice_distribution
9. �.�. �������� � �� "���������� ��������� � �������� � ���", �.,������������, 2004 %������������ ?????????? - �� ������������ (����)
10. �.�. ������, �.�. ��������� � ��. "�����, ����������� � ��������� ���������� �������� � ����������. ���. (���)", �., "���. �����", 1975, 296�. %������������ / �������� ��������
11. ���� �� ������� � �������� % ��� ����� ����������
12. ����� ������ ��� �������������� � fast � slow ������������
13. Lecture notes on state estimation of nonlinear non-Gaussian stochastic systems %���������� ������� �������-��� [12]
14. ���-������ �� ���������� �-��� ������� � ����
15. PSIAKI,2002, EKF methods for traking weak GPS signals
16.
%}



%}