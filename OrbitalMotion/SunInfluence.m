function [ influence ] = SunInfluence( time, coordinates )
    %% входные данные:
    %T_in - текущее время, на которое необходимо вычислить воздействия Солнца и Луны
    %Coord - массив координат. Coord(1) - x; Coord(2) - y; Coord(3) - z
    %Выходные данные
    %Sun - массив ускорений от Солнца, 1-x; 2-y; 3-z;
    %Moon - массив ускорений от Луны, 1-x; 2-y; 3-z;
    %%
    % Constant
    muSun = 0.1325263e12; %[км^3/с^2] - константа гравитационного поля Солнца
    
    aSun = 1.49598e8;%большая полуось «орбиты» Солнца
    
    eSun = 0.016719;%эксцентриситет солнечной «орбиты»
    
    qSun = 6.2400601269 + 628.3019551714*time - 0.0000026820*time^2;%[рад] - Средняя аномалия Солнца
    
    wSun = -7.6281824375 + 0.0300101976*time + 0.0000079741*time^2;%[рад] - Средняя тропическая долгота перигея орбиты Солнца,
    
    epsilon = 0.4090926006 - 0.0002270711*time;%[рад] - Средний наклон эклиптики к экватору
    
    %% Kepler's equation for Sun
    excentrAnomSun = qSun;
    i = 2;
    while abs( qSun+eSun*sin( excentrAnomSun ) - excentrAnomSun ) >= 1e-8
        excentrAnomSun = qSun+eSun*sin( excentrAnomSun );
        i=i+1;
    end;
    
    %% Other constants    
    sinFiSun=( sqrt(1-eSun^2)*sin(excentrAnomSun) )/( 1-eSun*cos(excentrAnomSun) );
    cosFiSun=( cos(excentrAnomSun) - eSun )/( 1-eSun*cos(excentrAnomSun) );
    
    eSun=cosFiSun*cos(wSun) - sinFiSun*sin(wSun);
    ettaSun=( sinFiSun*cos(wSun) + cosFiSun*sin(wSun) )*cos(epsilon);
    fSun=( sinFiSun*cos(wSun) + cosFiSun*sin(wSun) )*sin(epsilon);
    rSun=aSun*( 1-eSun*cos(excentrAnomSun) );
    
    %% Calculation of accelerations    
    xNormSun = coordinates( 1 ) / rSun;
    yNormSun = coordinates( 2 ) / rSun;
    zNormSun = coordinates( 3 ) / rSun;
    muNormSun = muSun / rSun^2;
    deltaSun = sqrt( ( eSun - xNormSun)^2 + ( ettaSun - yNormSun )^2 + ( fSun - zNormSun )^2 );
    
    influence(1, 1) = muNormSun*( ( eSun - xNormSun )/deltaSun^3 - eSun );
    influence(2, 1) = muNormSun*( ( ettaSun - yNormSun )/deltaSun^3 - ettaSun );
    influence(3, 1) = muNormSun*( ( fSun - zNormSun )/deltaSun^3 - fSun );
end
