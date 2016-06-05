%% функция расчёта ускорений от Солнца на заданное время закладки эфемерид
% см ИДК ГЛОНАСС П.3.1.1
%aster==asterisk==* - в именах переменных
function [Sun] = SunInfluence(T_in, Coord)
%% входные данные:
%T_in - текущее время, на которое необходимо вычислить воздействия Солнца и Луны
%Coord - массив координат. Coord(1) - x; Coord(2) - y; Coord(3) - z
%Выходные данные
%Sun - массив ускорений от Солнца, 1-x; 2-y; 3-z;
%Moon - массив ускорений от Луны, 1-x; 2-y; 3-z;
%% 
% Constant
    Mu_sun = 0.1325263e12; %[км^3/с^2] - константа гравитационного поля Солнца

    A_sun = 1.49598e8;%большая полуось «орбиты» Солнца

    e_sun = 0.016719;%эксцентриситет солнечной «орбиты»

    q_sun = 6.2400601269 + 628.3019551714*T_in - 0.0000026820*T_in^2;%[рад] - Средняя аномалия Солнца

    W_sun = -7.6281824375 + 0.0300101976*T_in + 0.0000079741*T_in^2;%[рад] - Средняя тропическая долгота перигея орбиты Солнца,

    Epsilon = 0.4090926006 - 0.0002270711*T_in;%[рад] - Средний наклон эклиптики к экватору

    %% Kepler's equation for Sun
    Excentr_anom_sun = q_sun;
    i = 2;
    while abs( q_sun+e_sun*sin( Excentr_anom_sun ) - Excentr_anom_sun ) >= 1e-8
        Excentr_anom_sun = q_sun+e_sun*sin( Excentr_anom_sun );    
        i=i+1;
    end;

    %% Other constants

    Sin_Fi_sun=( sqrt(1-e_sun^2)*sin(Excentr_anom_sun) )/( 1-e_sun*cos(Excentr_anom_sun) );
    Cos_Fi_sun=( cos(Excentr_anom_sun) - e_sun )/( 1-e_sun*cos(Excentr_anom_sun) );

    E_sun=Cos_Fi_sun*cos(W_sun) - Sin_Fi_sun*sin(W_sun);
    Etta_sun=( Sin_Fi_sun*cos(W_sun) + Cos_Fi_sun*sin(W_sun) )*cos(Epsilon);
    F_sun=( Sin_Fi_sun*cos(W_sun) + Cos_Fi_sun*sin(W_sun) )*sin(Epsilon);
    r_sun=A_sun*( 1-e_sun*cos(Excentr_anom_sun) );
    %-----------------------------------------------------------
    %% Calculation of accelerations

    X_norm_sun = Coord( 1 ) / r_sun;
    Y_norm_sun = Coord( 2 ) / r_sun;
    Z_norm_sun = Coord( 3 ) / r_sun;
    Mu_norm_sun = Mu_sun / r_sun^2;
    delta_sun = sqrt( ( E_sun - X_norm_sun)^2 + ( Etta_sun - Y_norm_sun )^2 + ( F_sun - Z_norm_sun )^2 );

    Sun(1,1) = Mu_norm_sun*( ( E_sun - X_norm_sun )/delta_sun^3 - E_sun );
    Sun(2,1) = Mu_norm_sun*( ( Etta_sun - Y_norm_sun )/delta_sun^3 - Etta_sun );
    Sun(3,1) = Mu_norm_sun*( ( F_sun - Z_norm_sun )/delta_sun^3 - F_sun );
end