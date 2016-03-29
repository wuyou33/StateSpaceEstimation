function [out] = Calc_param(satpos_xyz_Rec, satpos_xyz_gln, alm_gln, ns)
%Функция Calc_param производит вычисления геометрической и радиовидимости,
%задержки сигнала, вычисление доплеровского сдвига и геометрического
%фактора.
%Входные данные:
%satpos_xyz_Rec - структура данных вектора положения КА
%satpos_xyz_gln - структура данных вектора положения НКА
%alm_gln        - структура альманаха НКА%ns
%ns             - номер НС

%Выходные данные:
%out.vis - параметр определяющий относительную видимость
%out.tz  - время задержки
%out.fdop- доплеровский сдвиг
%out.cH  - косинусы

%========================Константы / Constants=============================
%{
%ДН антенн, используемых на КА и НКА

%ГЛОНАСС
DN_GLN = [0 10.7; 9 11.763636;10 11.881818; 11 12; 12 12; 15 12; 19 9;...
          20 7; 21 5.5; 22 4; 23 2.7; 24 1.5; 25 0; 30 -5; 35 -4.375;...
          40 -5.875; 45 -7.875; 50 -8.25; 55 -10.875; 60 -12.75];
%GPS - в расчетах не используется
DN_GPS = [0 13; 5 13.8; 6 14; 8 14.7; 9 14.9; 10 15; 11 14.9; 12 14.7;...
          14 13.7; 18 6; 20 0; 22 -10; 23 -8; 24 -5; 26 -3; 30 -2; 34 -3;...
          36 -5; 38 -10; 40 -15];
%КА
DN_ARN = [0 3; 2 3; 4 3; 6 3; 8 3; 10 3; 12 2.96; 14 2.94; 16 2.92; 18 2.935;...
          20 2.925; 22,2.9235; 24 2.914; 29 2.9; 60 1.5; 90 -2];
%}

%Частоты сигналов
f0(1) = 1602*10^6;          %частота 1 диапазона
f0(2) = 1246*10^6;          %частота 2 диапазона
df(1) = 562.5*10^3;         %дискрет 1 диапазона
df(2) = 437.5*10^3;         %дискрет 2 диапазона

%==========================================================================


fL(1) = f0(1) + alm_gln.Nn(ns) * df(1);
fL(2) = f0(2) + alm_gln.Nn(ns) * df(2);

[vision,P_input, cH] = Radio_vision_satelate(satpos_xyz_Rec, satpos_xyz_gln);

for i = 1:1%2 для L2 диапазона
   
    if (i == 1)
        [fdop1, tz, N_Chip1,N_Navig1,fi1, t_chip_resid1,N_Fdop_chip1,T_navig_resid] = signal_param(satpos_xyz_Rec, satpos_xyz_gln, fL(i));
%         fdop1 = fdop;
%         fi1   = fi;
        
    elseif (i == 2)
%         [fdop, tz, fi, t_chip_resid2] = signal_param(satpos_xyz_Rec, satpos_xyz_gln, fL(i));
%         fdop2 = fdop;
%         fi2   = fi;
    end
    
end

out.vis     = vision;
out.Power=P_input;
out.Tz      = tz;
out.N_navig=N_Navig1;
out.N_chip_Tau=N_Chip1;
out.T_chip_resid=t_chip_resid1;
out.T_navig_resid=T_navig_resid;
% out.T_chip_resid_2=t_chip_resid2;
out.Fdop_1  = fdop1;
% out.Fdop_2  = fdop2;
out.N_chip_Fdop=N_Fdop_chip1;

out.Phi_HF_1    = fi1;
% out.fi_2    = fi2;
out.cH      = cH;
out.Litera=alm_gln.Nn(ns);
out.delta_fi0=0;
out.d_fd_dt=0;
