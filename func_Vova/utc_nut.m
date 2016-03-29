function nut = utc_nut(t) 
%Имя: utc_nut 
%Функция предназначена для расчета нутации  
%Входные данные: 
%t =6.023472005475702e+002; 
%Выходные данные: 
%nut - нутация 
   R = 1296000; % ( 1r=360grad=1 296 000 cek)  
    RAD_SEK_ANGL   =  pi/(3600*180); 
    t2 = t * t; 
    t3 = t2 * t; 
    l = 485866.733 + (1325.0 * R + 715922.633) * t + 31.310 * t2 + 0.064 * t3;%1.034807679476340e+012 
    l1 = 1287099.804 + (99 * R + 1292581.224) * t - 0.577 * t2 - 0.012 * t3; 
    f = 335778.877 + (1342 * R + 295263.137) * t - 13.257 * t2 + 0.011 * t3; 
    dd = 1072261.307 + (1236 * R + 1105601.328) * t - 6.891 * t2 + 0.019 * t3; 
    omega = 450160.280 - (5 * R + 482890.539) * t + 7.455 * t2 + 0.008 * t3; 
    eps0 = 84381.448 - 46.8150 * t - 0.00059 * t2 + 0.001813 * t3; 
%    eps0 = 84381.447996 - 46.8150 * t - 0.00059 * t2 + 0.001813 * t3; 
    eps_d = utc_nut_fi_eps(t, l, l1, f, dd, omega, 'd','e');  
    eps_k = utc_nut_fi_eps(t, l, l1, f, dd, omega, 'k','e'); 
    eps = eps0 + eps_d + eps_k;  
    cos_eps = cos(RAD_SEK_ANGL * eps ) / 15.0; 
    d_fi = utc_nut_fi_eps(t, l, l1, f, dd, omega, 'd', 'f'); 
    k_fi = utc_nut_fi_eps(t, l, l1, f, dd, omega, 'k', 'f'); 
    nut1 = d_fi * cos_eps; 
    nut2 = k_fi * cos_eps; 
% nut3 = 0.00264 * sin(omega) + 0.000063 * sin(2.0 * omega)  
   nut3 = 0; 
   nut1_dop = nut1; 
   nut2_dop = nut2; 
   nut3_dop = nut3; 
   nut = nut1 + nut2 + nut3; 
