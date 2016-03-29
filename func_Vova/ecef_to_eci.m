function [satpos_eci] = ecef_to_eci(S, satpos_ecef) 
%»м€ функции:eci_to_ecef 
%‘ункци€ преобразовани€ координат 
%¬ходные данные: 
%S - угловое положение неинерциальной —  относительно инерциальной
%—труктура satpos_eci 
%satpos_ecef.x -  координата x в подвижной системе координат (ECEF); 
%satpos_ecef.y - координата  y в подвижной системе координат (ECEF); 
%satpos_ecef.z - координата  z в подвижной системе координат (ECEF); 
%satpos_ecef.vx - скорость vx в подвижной системе координат (ECEF); 
%satpos_ecef.vy - скорость vy в подвижной системе координат (ECEF); 
% satpos_ecef.vz - скорость vz в подвижной системе координат (ECEF);

%¬ыходные данные: 
% —труктура satpos_eci  
%satpos_eci.x - координата x в абсолютной неподвижной системе координат (ECI); 
%satpos_eci.y - координата y в абсолютной неподвижной системе координат (ECI);
%satpos_eci.z - координата z в абсолютной неподвижной системе координат (ECI);
%satpos_eci.vx - скорость по оси x в абсолютной неподвижной системе координат (ECI); 
%satpos_eci.vy - скорость по оси z в абсолютной неподвижной системе координат (ECI);
%satpos_eci.vz-  скорость по оси z в абсолютной неподвижной системе координат (ECI); 
% оэффициенты 
global GL_W_rot_Earth
    A_trans=[cos(S) -sin(S) 0;...
                    sin(S)  cos(S) 0;...
                       0          0     1];% матрица поворота
              
    W=[0 -GL_W_rot_Earth 0;GL_W_rot_Earth 0 0;0 0 0];%матрица угловых скоростей
    
    
    oxyz = [satpos_ecef.x; satpos_ecef.y; satpos_ecef.z];
    voxyz = [satpos_ecef.vx; satpos_ecef.vy; satpos_ecef.vz];
    
    
    outxyz=A_trans*oxyz;
    voutxyz=A_trans*voxyz+W*outxyz;
       
    satpos_eci.x =  outxyz(1); 
    satpos_eci.y =  outxyz(2); 
    satpos_eci.z =  outxyz(3);
    %доработать со скорост€ми
    satpos_eci.vx = voutxyz(1);
    satpos_eci.vy = voutxyz(2);
    satpos_eci.vz = voutxyz(3);
