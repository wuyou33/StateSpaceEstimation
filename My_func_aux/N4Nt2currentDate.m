%функция перевода данных о дате в навигационном сообщении в формат dd.mm.yyyy
%алгоритм взят из ИДК Глонасс П.3.1.3
function [Day, Month, Year]=N4Nt2currentDate(n4,nt)

JD =1461*(n4-1) + nt + 2450082.5;
JDN = JD+ 0.5;
a = JDN+ 32044;
b = fix( (4*a + 3) / 146097 );
c = a - fix( (146097*b) / 4 );
d = fix( (4*c + 3) / 1461 );
e = c - fix( (1461*d) / 4 );
m = fix( (5*e + 2) / 153 );

Day = e - fix( (153*m + 2) / 5 ) + 1;
Month = m + 3 - 12*fix(m / 10);
Year = 100*b + d - 4800 + fix(m / 10);