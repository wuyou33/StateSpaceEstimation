function [ND]=SNC(n_i,NS,Num,N0,ti,P3,alm_gln,coord_gln);
%Функция заполнения Кадра, навигационными данными
%Входные данные:
%n_i - номер временного интервала
%NS  - номер спутника
%Num - номер кадра
%N0  - номер суток создания кадра
%ti  - текущее время
%P3  - признак кадра(0 - обычный(1-4й кадры), 1 - 5й кард)
%Выходные данные
%
init_alm = 1 + ((Num-1)*5);     %начальное значение альманаха спутников

    %1-я строка
    ND(1).m  = 1;
        %2 бита пустые
    ND(1).P1 = 01;%признак смены ОИ (const = 30мин)
    ND(1).tk = tk_calc(ti);
    ND(1).vx = coord_gln(NS,n_i).vx;
    ND(1).ax = 1;
    ND(1).x  = coord_gln(NS,n_i).x;
    
    %2-я строка
    ND(2).m  = 2;
    ND(2).Bn = 000; %const - что спутник пригоден для связи
    ND(2).P2 = 1;   %const - tb - нечетная
    ND(2).tb = floor((ti-900)/1800)+1;%алгоритм расчета для этой цифры
        %пропуск 5 бит
    ND(2).vy = coord_gln(NS,n_i).vy;
    ND(2).ay = 0;
    ND(2).y  = coord_gln(NS,n_i).y;
    
    %3-я строка
    ND(3).m  = 3;
    ND(3).P3 = P3;%1 or 0 альманах 5 или 4 нка
    ND(4).yntb= 0.001;%относительно отклонение значения несущей
        %пропуск 1 бит
    ND(3).p  = 1;%предстоит разобраться
    ND(3).ln = alm_gln.Cn(NS);%пригодность нка
    ND(3).vz = coord_gln(NS,n_i).vz;
    ND(3).az = 0;
    ND(3).z  = coord_gln(NS,n_i).z;
    
    %4-я строка
    ND(4).m   = 4;
    ND(4).tntb= alm_gln.Tgpgl(NS);%сдвиг шкалы времени нка относительно шкалы времени глонасс
    ND(4).dtn = 0;%смещение сигнала tf2-tf1
    ND(4).En  = N0 - alm_gln.Na(NS);
        %пропуск 14 бит
    ND(4).P4  = 0;%надо разобраться
    ND(4).Ft  = 0;%фактор точности измерений
        %пропуск 3 бита
    ND(4).Nt  = alm_gln.Na(NS);%текущая дата в 4х летнем интервале
    ND(4).n   = alm_gln.n(NS);%номер НКА
    ND(4).M   = 01;%модификация спутника - Глонасс-М
    
    %5-я строка
    ND(5).m   = 5;
    ND(5).Na  = alm_gln.Na(NS);
    ND(5).tc  = alm_gln.Tglutc(NS);
        %1 bit
    ND(5).N4  = floor((alm_gln.date.year-1996)/4);%надо вычислять но сейчас он равен 6
    ND(5).tgps= alm_gln.Tgpgl(NS);
    ND(5).ln  = alm_gln.Cn(NS);%пригодность нка
    
%%%%%%%%%%%%%%%%%%%========================================Альманах системы

    %6-я строка - альманах спутников
    ND(6).m   = 6;
    ND(6).Cn  = alm_gln.Cn(init_alm);
    ND(6).Man = alm_gln;%признак модификации
    ND(6).na  = alm_gln.n(init_alm);
    ND(6).tan = alm_gln.Tka(init_alm);
    ND(6).lam = alm_gln.Lam(init_alm);
    ND(6).di  = alm_gln.dI(init_alm);
    ND(6).ecc = alm_gln.E(init_alm);

    %7-я строка
    ND(7).m   = 7;
    ND(7).wn  = alm_gln.omegan(init_alm);
    ND(7).tlam= alm_gln.tc(init_alm);
    ND(7).dT  = alm_gln.dT(init_alm);
    ND(7).dTT = alm_gln.dTT(init_alm);
    ND(7).Ha  = alm_gln.Na(init_alm);
    ND(7).ln  = alm_gln.Cn(init_alm);
    
    %8-я строка альманах спутника init_alm+1
    ND(8).m   = 8;
    ND(8).Cn  = alm_gln.Cn(init_alm+1);
    ND(8).Man = alm_gln;%признак модификации
    ND(8).na  = alm_gln.n(init_alm+1);
    ND(8).tan = alm_gln.Tka(init_alm+1);
    ND(8).lam = alm_gln.Lam(init_alm+1);
    ND(8).di  = alm_gln.dI(init_alm+1);
    ND(8).ecc = alm_gln.E(init_alm+1);

    %9-я строка
    ND(9).m   = 9;
    ND(9).wn  = alm_gln.omegan(init_alm+1);
    ND(9).tlam= alm_gln.tc(init_alm+1);
    ND(9).dT  = alm_gln.dT(init_alm+1);
    ND(9).dTT = alm_gln.dTT(init_alm+1);
    ND(9).Ha  = alm_gln.Na(init_alm+1);
    ND(9).ln  = alm_gln.Cn(init_alm+1);
    
    %10-я строка альманах спутника init_alm+2
    ND(10).m   = 10;
    ND(10).Cn  = alm_gln.Cn(init_alm+2);
    ND(10).Man = alm_gln;%признак модификации
    ND(10).na  = alm_gln.n(init_alm+2);
    ND(10).tan = alm_gln.Tka(init_alm+2);
    ND(10).lam = alm_gln.Lam(init_alm+2);
    ND(10).di  = alm_gln.dI(init_alm+2);
    ND(10).ecc = alm_gln.E(init_alm+2);

    %11-я строка
    ND(11).m   = 11;
    ND(11).wn  = alm_gln.omegan(init_alm+2);
    ND(11).tlam= alm_gln.tc(init_alm+2);
    ND(11).dT  = alm_gln.dT(init_alm+2);
    ND(11).dTT = alm_gln.dTT(init_alm+2);
    ND(11).Ha  = alm_gln.Na(init_alm+2);
    ND(11).ln  = alm_gln.Cn(init_alm+2);
    
    %12-я строка альманах спутника init_alm+3
    ND(12).m   = 12;
    ND(12).Cn  = alm_gln.Cn(init_alm+3);
    ND(12).Man = alm_gln;%признак модификации
    ND(12).na  = alm_gln.n(init_alm+3);
    ND(12).tan = alm_gln.Tka(init_alm+3);
    ND(12).lam = alm_gln.Lam(init_alm+3);
    ND(12).di  = alm_gln.dI(init_alm+3);
    ND(12).ecc = alm_gln.E(init_alm+3);

    %13-я строка
    ND(13).m   = 13;
    ND(13).wn  = alm_gln.omegan(init_alm+3);
    ND(13).tlam= alm_gln.tc(init_alm+3);
    ND(13).dT  = alm_gln.dT(init_alm+3);
    ND(13).dTT = alm_gln.dTT(init_alm+3);
    ND(13).Ha  = alm_gln.Na(init_alm+3);
    ND(13).ln  = alm_gln.Cn(init_alm+3);
    
    %14-я строка альманах спутника init_alm+4
    ND(14).m   = 14;
    ND(14).Cn  = alm_gln.Cn(init_alm+4);
    ND(14).Man = alm_gln;%признак модификации
    ND(14).na  = alm_gln.n(init_alm+4);
    ND(14).tan = alm_gln.Tka(init_alm+4);
    ND(14).lam = alm_gln.Lam(init_alm+4);
    ND(14).di  = alm_gln.dI(init_alm+4);
    ND(14).ecc = alm_gln.E(init_alm+4);

    %15-я строка
    ND(15).m   = 15;
    ND(15).wn  = alm_gln.omegan(init_alm+4);
    ND(15).tlam= alm_gln.tc(init_alm+4);
    ND(15).dT  = alm_gln.dT(init_alm+4);
    ND(15).dTT = alm_gln.dTT(init_alm+4);
    ND(15).Ha  = alm_gln.Na(init_alm+4);
    ND(15).ln  = alm_gln.Cn(init_alm+4);