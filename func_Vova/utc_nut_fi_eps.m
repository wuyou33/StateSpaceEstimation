 function nut_fi_eps = utc_nut_fi_eps(t, l, l1, f, dd, omega, typ_nut, fi_eps) 
% Имя функции:utc_nut_fi_eps 
%Функция предназначена для расчета нутации по епсилон и  фи 
%Входные данные: 
%t - значение параметра приведено в функции utc_nut; 
%l - значение параметра приведено в функции utc_nut; 
%l1 - значение параметра приведено в функции utc_nut; 
%f -значение параметра приведено в функции utc_nut; 
%dd - значение параметра приведено в функции utc_nut; 
% оmega - значение параметра приведено в функции utc_nut; 
%typ_nut - признак параметра приведен в функции utc_nut; 
%fi_eps - признак параметра приведен в функции utc_nut; 
%Выходные данные: 
%nut_fi_eps - значение нутации по по епсилон и  фи 
% применить функцию koef 
[koef_id, koef_abd, koef_ik, koef_abk] = koef; 
 RAD_SEK_ANGL  =  pi/(3600*180); 
    if (typ_nut == 'd') 
        n = 30; 
    else 
          n = 76; 
    end; 
    sum_a = 0; 
    sum_b = 0; 
   for i = 1 : n 
   if (typ_nut == 'd')  
   s1 = koef_id(i,1) * l + koef_id(i,2) * l1 + koef_id(i,3) * f + koef_id(i,4) * dd + koef_id(i,5) * omega; 
        if (fi_eps == 'f')  
             a  = koef_abd(i,1) * 1e-4; 
           bt = koef_abd(i,2) * 1e-4; 
            else 
            a  = koef_abd(i,3) * 1e-4; 
           bt = koef_abd(i,4) * 1e-4; 
            end; 
    else  
    s1 = koef_ik(i,1) * l + koef_ik(i,2) * l1 + koef_ik(i,3) * f + koef_ik(i,4) * dd + koef_ik(i,5) * omega; 
     if (fi_eps == 'f')  
             a  = koef_abk(i,1) * 1e-4; 
           bt = koef_abk(i,2) * 1e-4; 
            else 
            a  = koef_abk(i,3) * 1e-4; 
           bt = koef_abk(i,4) * 1e-4; 
            end; 
    end; 
      if (fi_eps == 'f')  
     sin_s1 = sin(RAD_SEK_ANGL * s1); 
        sa = a * sin_s1; 
     sb = bt * sin_s1; 
    else  
       cos_s1 = cos(RAD_SEK_ANGL * s1); 
       sa = a  * cos_s1; 
       sb = bt * cos_s1; 
     end; 
     arg = RAD_SEK_ANGL * s1; 
  sum_a = sum_a + sa; 
  sum_b = sum_b + sb; 
end; 
nut_fi_eps = sum_a + sum_b * t; 
