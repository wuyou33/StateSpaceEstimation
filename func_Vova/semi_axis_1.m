% double   semi_axis(double t_dr, double i_incl, double  ecc, double  omega_n) 
function [a_npp] = semi_axis_1(t_dr,  i_incl,   ecc,   omega_n) 
%Имя функции:semi_axis_1 
%Функция вычисляет радиус орбиты спутника ГЛОНАСС в соответствии с интерфейсным контроль-ным  
%документом ГЛОНАСС 
%Входные данные: 
%t_dr - драконический период обращения спутника ГЛОНАСС (секунды) 
%i_incl - наклонение орбиты спутника ГЛОНАСС (радиан)  
%ecc - эксцентриситет 
%omega_n - аргумент перигея орбиты спутника ГЛОНАСС (радиан) 
%Выходные данные: 
%a_npp - большая полуось орбиты спутника ГЛОНАСС (километр) 
global GL_Ae %[км] - экваториальный радиус Земли
global J2_0;% - вторая зональная гармоника геопотенциала
global GL_MuE%[км^3/с^2] - константа гравитационного поля Земли

% B_PZ90_KM = 6356.75136174571344; % AP_LAND (Km) Polar radius of the Earth */ 

    epsilon = 10.0e-3; 
    sin_i = sin(i_incl); 
    sin_i2 = sin_i * sin_i; 
   % v = -omega_n;% omega_n - angle of a perigee %?????????
    ecc2_1 = 1.0 - ecc * ecc; 
    b1 = 2.0 - 2.5 * sin_i2; 
    b2 = sqrt(ecc2_1 * ecc2_1 * ecc2_1 ); 
    b3 = 1.0 + ecc * cos(omega_n*pi); %?????????????????
    b4 = b3 * b3 * b3;
    b3 = b3 * b3; 
    b5 = b4 / ecc2_1; 
    b = b1 * b2 / b3 + b5;      
    tock_2pi = t_dr./(pi * 2); 
    a_dop = GL_MuE.* tock_2pi.* tock_2pi; 
    a_n = a_dop.^(1/3); 

%     nn = 0; 
 a_npp=zeros(1,length(t_dr));
for j=1:length(t_dr)
 
 dda = epsilon + 1; 
 while ( dda > epsilon ) 
       p = a_n(j) * ecc2_1; % Focal parameter   
       ae_p = GL_Ae / p; 
       b0 = 1.0 - (1.5 * J2_0 * ae_p * ae_p) * b; %???
       t_ock = t_dr(j) / b0; 
       tock_2pi = t_ock / (pi * 2); 
       a_dop2 = GL_MuE * tock_2pi * tock_2pi; 
       a_npp(j) = a_dop2^(1/3); 
       dda = abs(a_n(j) - a_npp(j)); 
       a_n(j) = a_npp(j); 
%        nn = nn + 1; 
 end; 
 
end;
 