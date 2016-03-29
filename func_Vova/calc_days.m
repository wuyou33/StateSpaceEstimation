function [days] = calc_days(DnMon , date , days_year);
%Под функция вычисления текущего номера суток.
%Вход:
%DnMon     - массив с кол-вом дней в месяцах
%date      - номер дня
%days_year - кол-во дней прошедших с последнего високосного года
%Выход:
%days      - номер текущих суток от последнего високосного года

days = 0;
%days_year = 0;

for i=1:(date.mon-1)
    
    days = days+DnMon(i);

end

days =days_year + days + date.day;