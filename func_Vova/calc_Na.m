function [days] = calc_Na(date)
%‘нкци€ вычислени€ календарного номера суток из текущей даты.
%¬ходные данные:
%—труктура date - в формате день.мес€ц.год
%date.day    -день
%date.mon    -мес€ц
%date.year   -год
%¬ыходные данные:
%Na -календарный номер суток от ближайшего високосного года
% clc;
% clear;
% alm_gln.date.day  = 10;       %тестовые данные
% alm_gln.date.mon  = 01;
% alm_gln.date.year = 2014;


DnMonV = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
DnMon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 

year_1 = mod(date.year, 4);
leap_year = date.year - year_1;

switch (year_1)
    case (1)
    days_year = 366;
    case (2)
    days_year = 731;
    case (3)
    days_year = 1096;
    otherwise
    days_year = 0;
end

if (year_1 == 0)
    days = calc_days(DnMonV , date , days_year);
else
    days = calc_days(DnMon , date , days_year);
end