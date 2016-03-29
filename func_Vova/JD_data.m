function [JD, day_year] = JD_data(timeUTC) 
%Имя:JD_data 
% Функция JD_data(timeUTC)  вычисляет : 
%JD - номер юлианского дня, day_year - номер дня года. 
%Входные данные: 
%Структура timeUTC  
%timeUTC.year - год, 
% timeUTC.mon - месяц, 
% timeUTC.day - день. 
%Выходные данные:  
%JD - юлианский день; 
%day_year- день от начала года. 
%количество дней в месяцах 
  DnMon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
%Вычисление номера юлианского дня опорной эпохи 
  jd0 = JD_epohi(timeUTC.year); 
%Учет високосного года 
nfebr = 0; 
if mod(timeUTC.year,4) == 0 
    nfebr = 1; 
end; 
%Расчетномера дня года 
   k = 0; 
   for i = 2 : timeUTC.mon 
        k = k + DnMon(i - 1); 
      if (i == 2)  
            k = k + nfebr; 
        end;   
   end; 
    day_year = k + timeUTC.day; 
%Расчет номера юлианского дня 
    JD = jd0 + day_year; 
