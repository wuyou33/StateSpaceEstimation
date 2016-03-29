%function [rin]=pz2oxyz(PathToTheFileRine);
%
%Функция пересчета координат
%из системы ПЗ-90-02 Oxyz в абсолютную Ox0y0z0
%
%
%
%
%
%
timeUTC.day = 6; 
timeUTC.mon = 09;
timeUTC.year= 2001;
nut     = 1;            %Признак вычисления истенного звездного времени
te=1;
S = s0_Nut(timeUTC, nut);

%Преобразование координат
x0(te) = x(te)*cos(S(te)) - y(te)*sin(S(te));
y0(te) = x(te)*sin(S(te)) + y(te)*cos(S(te));
z0(te) = z(te);

%Преобразование скоростей
Vx0(te) = Vx(te)*cos(S(te)) - Vy(te)*sin(S(te)) - Wz*y0(te);
Vy0(te) = Vx(te)*cos(S(te)) + Vy(te)*cos(S(te)) + Wz*x0(te);
Vz0(te) = Vz(te);
