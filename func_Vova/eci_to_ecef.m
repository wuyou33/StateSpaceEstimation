function [satpos_ecef] =eci_to_ecef(s0, ti, satpos_eci) 
%Имя функции:eci_to_ecef 
%Функция преобразования координат 
%Входные данные: s0 - истинное звездное время в текущий момент обсервации , 
%ti - текущее время; satpos_eci  
%Структура satpos_eci 
%satpos_eci.x -  координата x в абсолютной неподвижной системе координат (ECI); 
%satpos_eci.y - координата  y в абсолютной неподвижной системе координат (ECI); 
%satpos_eci.z - координата  z в абсолютной неподвижной системе координат (ECI); 
%satpos_eci.vx - скорость vx в абсолютной неподвижной системе координат (ECI);  
%satpos_eci.vy - скорость vy в абсолютной неподвижной системе координат (ECI); 
% satpos_eci.vz - скорость vz в абсолютной неподвижной системе координат (ECI); 
%Выходные данные: 
% Структура satpos_ecef  
%satpos_ecef.x - координата x в подвижной системе координат (ECEF); 
%satpos_ecef.y - координата y в подвижной системе координат (ECEF); 
%satpos_ecef.z - координата z в подвижной системе координат (ECEF); 
%satpos_ecef.vx - скорость по оси x в подвижной системе координат (ECEF); 
%satpos_ecef.vy - скорость по оси z в подвижной системе координат (ECEF); 
%satpos_ecef.vz-  скорость по оси z в подвижной системе координат (ECEF); 
%Коэффициенты 
% SEC_IN_RAD - коэффициент преобразования секунд в радианы 
%  s0(radian) = s0 (sek) * SEC_IN_RAD, where 
 % SEC_IN_RAD = 2 * pi / (24 * 3600)  = pi / 43200 
 SEC_IN_RAD = pi / 43200; 
 OMEGA_Z  = 0.7292115e-4;  %( скорость  вращения  Земли (angular speed of rotation of the Earth, рад/cek)  
    s_zv = s0 * SEC_IN_RAD + OMEGA_Z * ti; 
    cos_s = cos(s_zv); 
    sin_s = sin(s_zv); 
 
    satpos_ecef.x =  satpos_eci.x * cos_s + satpos_eci.y * sin_s; 
    satpos_ecef.y = -satpos_eci.x * sin_s + satpos_eci.y * cos_s; 
    satpos_ecef.z =  satpos_eci.z; 
 
    satpos_ecef.vx =  satpos_eci.vx * cos_s + satpos_eci.vy * sin_s + OMEGA_Z * satpos_ecef.y; 
    satpos_ecef.vy = -satpos_eci.vx * sin_s + satpos_eci.vy * cos_s - OMEGA_Z * satpos_ecef.x; 
    satpos_ecef.vz =  satpos_eci.vz; 
