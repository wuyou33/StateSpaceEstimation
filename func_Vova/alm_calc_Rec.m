function [sat_pos_xyz_Rec sat_vel_xyz_Rec] = alm_calc_Rec(tau, alm_test)
%Функция расчета вычисления координат приемника, движущегося по
%элиптической орбите
%Входные данные: альманах системы
%Выхдные данные: координаты и скорости на заданное время расчета
global GL_MuE%[км^3/с^2] - константа гравитационного поля Земли
   
i =(alm_test.I/180)*pi;%[рад] - угол наклонения орбиты
ecc_2 = alm_test.E*alm_test.E;%эксцентриситет в квадрате
% cos_i2 = cos(i)*cos(i);
% sin_i2 = sin(i)*sin(i);
   
%Вычисление относительного времени движения КА
dtpr = (tau - alm_test.tc);
%Вычисление фокального параметра
p = alm_test.a * (1 - ecc_2);
%Вычисление лямды
Lamda = sqrt(GL_MuE)/(alm_test.a^(3/2));
%Вычисление средней аномалии
M = Lamda*(dtpr) ;

E_pp = M;
E_npp = M + alm_test.E*sin(E_pp);
%Решение уравнения Кеплера - нахождение эксцентрической аномалии
while (abs( E_npp - E_pp )>=10^-9)
    E_pp = E_npp;    
    E_npp = M + alm_test.E*sin(E_pp);
end;
%Нахождение истенной аномалии произвольной точки
v = 2*atan2( sqrt(1+alm_test.E)*tan(E_npp/2),sqrt(1-alm_test.E));
%Нахождение углового расстояния от восходящего узла
u = v + alm_test.omegan;                            

r = p/(1+alm_test.E*cos(v));

v_r = sqrt(GL_MuE/p)*alm_test.E*sin(v);   %alm_test.omegan !!!!!!!!!!!!!!
v_u = sqrt(GL_MuE/p)*( 1+alm_test.E*cos(v) );%alm_test.omegan !!!!!!!!!!

m_x = cos(alm_test.Lam)*cos(u)-sin(alm_test.Lam)*sin(u)*cos(i);
m_y = sin(alm_test.Lam)*cos(u)+cos(alm_test.Lam)*sin(u)*cos(i);
m_z = sin(u)*sin(i);

sat_pos_xyz_Rec(1,1) = r*m_x;
sat_pos_xyz_Rec(2,1) = r*m_y;
sat_pos_xyz_Rec(3,1) = r*m_z;

sat_vel_xyz_Rec(1,1) = v_r*(cos(alm_test.Lam)*cos(u)-sin(alm_test.Lam)*sin(u)*cos(i))-...
                    v_u*(cos(alm_test.Lam)*sin(u)+sin(alm_test.Lam)*cos(u)*cos(i)); 

sat_vel_xyz_Rec(2,1) = v_r*(sin(alm_test.Lam)*cos(u)+cos(alm_test.Lam)*sin(u)*cos(i))-...
                    v_u*(sin(alm_test.Lam)*sin(u)-cos(alm_test.Lam)*cos(u)*cos(i)); 
    
sat_vel_xyz_Rec(3,1) = v_r*sin(u)*sin(i)+v_u*cos(u)*sin(i); 