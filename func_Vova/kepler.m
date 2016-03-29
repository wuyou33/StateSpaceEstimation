function [Ek] =  kepler(Mk, ecc, eps);   
%Имя функции:kepler 
%Функция предназначена для решения уравнения Кеплера 
%Вход:
%Mk -аномалия
%ecc-ексцесинтриситет
%eps-точность
%Выход:
%Ek -истенная аномалия 
%   eps = 1.0E-15; 
   y = ecc * sin(Mk); 
   x1 = Mk; 
   x = y; 
   for k = 0 : 15 
       x2 = x1; 
       x1 = x; 
       y1 = y; 
       y = Mk - (x - ecc * sin(x)); 
       if (abs(y - y1) < eps) 
           break 
           end  
      x = (x2 * y - x * y1) / (y - y1); 
    end  %k   
    Ek = x; 
