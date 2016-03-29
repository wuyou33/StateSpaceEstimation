function [fdop, tz, N_int_Chip,N_int_Navig,fi, t_chip_resid,N_int_Fdop_cells,T_navig_resid] = signal_param(rec, trans, fs)    
%Функция вычисляет дорплероский сдвиг, время задержки и фазу сигнала
%Вход:
%rec - структура координат КА
%tranc - структура координат НКА
%fs - значение частоты для данного НКА
%Выход:
%fdop - доплеровский свдиг частоты
%tz - время задержки
%fi - фаза сигнала

global GL_T_chip GL_T_Navig_message GL_Fdop_cell

N_total=length(rec.x);

D     = zeros(1,N_total);      %обнуление
S1    =  zeros(1,N_total);
S2    =  zeros(1,N_total);
c      = 299792.458;%[км/c] - скорость света в КМ


KA     = [rec.x;rec.y;rec.z];    
NKA    = [trans.x;trans.y;trans.z];

vKA    = [rec.vx; rec.vy; rec.vz];
vNKA   = [trans.vx; trans.vy; trans.vz];

%вычисление суммарного вектора
X12=zeros(3,N_total);
X21=zeros(3,N_total);

for i = 1:3
    DD = KA(i,:) - NKA(i,:);
    D  = D + DD.*DD;
    X12(i,:) = DD;            %вектор с направлением от НКА к КА
    X21(i,:) = -DD;           %вектор с направлением от КА к НКА
end;

D = sqrt(D);                %расстояние между объектами

    summ1 = sqrt(vKA(1,:).^2 + vKA(2,:).^2 + vKA(3,:).^2);     %длинна вектора скорости КА
    summ2 = sqrt(vNKA(1,:).^2 + vNKA(2,:).^2 + vNKA(3,:).^2);  %длинна вектора скорости НКА
    
    for i = 1:3
        S1 = S1+X12(i,:).*vNKA(i,:);  %скалярное произведение вектора визирования и вектора скорости приемника
    end;
    arg1 = acos(S1./(D.*summ2));  %вычисление угла поворота
    v1   = cos(arg1).*summ2;     %поворот вектора на нужный угол
    
    for i = 1:3
        S2 = S2+X21(i,:).*vKA(i,:); %скалярное произведение вектора визирования и вектора скорости передатчика
    end;
    arg2 = acos(S2./(D.*summ1));  %вычисление угла поворота
    v2   = cos(arg2).*summ1;     %поворот вектора на нужный угол

fdop = fs.*(v2 + v1)./c;          %Доплеровский сдвиг частоты

%=============================Задержка сигнала=============================
clear v2 v1 arg1 arg 2 S1 S2 summ1 summ2 X12 X21 KA NKA vKA vNKA

tz = D./c;

clear D
%==========================Фаза высокочастотного сигнала=====================
%{
Tp = 1/fs;

fi.int   = fix(tz/Tp);          %число целых периодов
t_resid       = mod(tz,Tp);          %оставшееся время в интервале 2 Pi
fi.res   = t_resid*2*180/Tp;               %то же, в градусах
%}    
% T_hf=1/fs; %период ВЧ
%==========================Фаза навигационного кода=====================

N_int_Navig = fix(tz/GL_T_Navig_message); %число целых периодов навигационного сообщения
N_int_Chip   = fix( (tz-N_int_Navig*GL_T_Navig_message)/GL_T_chip);          %число целых чипов дальномерного кода в нужном нам интервале по 1 мс!!!
T_navig_resid=mod(tz,GL_T_Navig_message);%время задержки от последнего целого навигациооного бита (<1 мс)
N_int_Fdop_cells=fix(fdop/GL_Fdop_cell);%число целых ячеек по 500 Гц Fдоп
t_chip_resid   = mod(tz,GL_T_chip);          %[сек] - оставшееся время после отбрасывания целого числа ЧИПОВ 
fi   = tz.*(fs+fdop)*2*pi;   %[рад] - полная фаза навигационного сообщения


