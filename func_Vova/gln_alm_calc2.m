function [calc] = gln_alm_calc2(ns,  n0, tmdv, time_s0, alm_gln);          
%Функция вычисления орбит спутников по ИКД для диапазона L3

MU = 398600.4418;
imid = 63;                  %среднее значение наклонения орбиты
Rz = 6378.136;
WZ = 7.2921150*10^-5;

%интервал времени с момента прохождения НКА 1го восходящего узла
deltat = (n0 - alm_gln.Na(ns)) *  86400 + tmdv - alm_gln.tc(ns);
%наклонение орбиты
i = (imid/180 + alm_gln.dI(ns))*pi;
%скорректированное движение
n = 2*pi / (43200+alm_gln.dT(ns)+alm_gln.dTT(ns)*deltat/(43200+alm_gln.dT(ns)));
%вычисление размера большой полуоси
a = (MU/(n*n))^(1/3);
%вычисление скорости изменения гринвичской долготы восходящего узла
OMEGAv = -10*sqrt((Rz/a)^7)*cos(i)/86400/180;
%аргумент перигея орбиты НКА
w = 5*sqrt((Rz/a)^7)*(5*(cos(i)*cos(i))-1)/86400/180;
%эксцентрическая аномалия
se = (sin(-( alm_gln.omegan(ns)+w*deltat)*pi))*sqrt(1-alm_gln.E(ns)*alm_gln.E(ns));
ce = alm_gln.E(ns) + cos(-( alm_gln.omegan(ns) +w*deltat)*pi);

Ew=atan2(se,ce);
% if ( ce > 0)
%     Ew = atan2(se,ce);
% elseif ((se>=0) && (ce<0))
%     Ew = pi+atan2(se,ce);
% elseif ((se<0) && (ce<0))
%     Ew = -pi+atan2(se,ce);
% elseif ((se>0) && (ce == 0))
%     Ew = pi/2;
% elseif ((se<0) && (ce == 0))
%     Ew = -pi/2;
% elseif ((se==0) && (ce == 0))
%     Ew = 0;
% end;

%вычисление средней аномалии на данный момент
deltati = (Ew - alm_gln.E(ns)*sin(Ew))/n;
E_0 = n*(deltat+deltati);

E_pp = E_0;
E_npp = E_pp+1;
   
while (abs( E_npp - E_pp )>=10^-9)
    E_pp = E_npp;    
    E_npp = E_0 - alm_gln.E(ns)*sin(Ew);
end;
E = E_npp;
sSig = sqrt(1 - alm_gln.E(ns)*alm_gln.E(ns))*sin(E);
cSig = cos(E) - alm_gln.E(ns);
Sig  = atan2(sSig,cSig);
%аргумент широты 
u = Sig + (alm_gln.omegan(ns)+w*deltat)*pi;
%модуль радиус вектора
r  = a*(1 - alm_gln.E(ns)*cos(E));
OMEGA = alm_gln.Lam(ns)*pi + (OMEGAv-WZ)*deltat;

ax = cos(u)*cos(OMEGA)-sin(u)*sin(OMEGA)*cos(i);
ay = cos(u)*sin(OMEGA)+sin(u)*cos(OMEGA)*cos(i);
az = sin(u)*sin(i);

calc.x = r*ax;
calc.y = r*ay;
calc.z = r*az;

calc.vx = 1;
calc.vy = 1;
calc.vz = 1;

   
    