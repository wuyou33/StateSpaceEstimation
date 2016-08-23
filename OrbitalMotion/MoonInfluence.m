function [ influence ] = MoonInfluence( time, coordinates )
    %
    % Constant
    muMoon = 4902.835; %[км^3/с^2] - константа гравитационного пол€ Ћуны
    
    aMoon = 3.84385243e5;%больша€ полуось орбиты Ћуны
    
    eMoon = 0.054900489;%эксцентриситет лунной орбиты,
    
    iMoon = 0.0898041080;%[рад] - —реднее наклонение орбиты Ћуны к плоскости эклиптики
    
    qMoon = 2.3555557435 + 8328.6914257190*time + 0.0001545547*time^2;%[рад] - —редн€€ аномали€ Ћуны
    
    omegaMoon = 2.1824391966 - 33.7570459536*time + 0.0000362262*time^2;%[рад] - —редн€€ долгота восход€щего узла Ћуны
    
    gMoon = 1.4547885346 + 71.0176852437*time - 0.0001801481*time^2;%[рад] - —редн€€ долгота периге€ орбиты Ћуны
    
    wpsilon = 0.4090926006 - 0.0002270711*time;%[рад] - —редний наклон эклиптики к экватору
    
    cosOmegaMoonOnSiniMoon = cos(omegaMoon)*sin(iMoon);
    ettaAster = sin(omegaMoon)*sin(iMoon);
    eAster = cos(iMoon);
    
    e11 = sin(omegaMoon)*cos(omegaMoon)*(1-cos(iMoon));
    e12 = 1-( sin(omegaMoon) )^2*(1-cos(iMoon));
    
    etta11=eAster*cos(wpsilon) - cosOmegaMoonOnSiniMoon*sin(wpsilon);
    etta12=e11*cos(wpsilon) + ettaAster*sin(wpsilon);
    
    f11=eAster*sin(wpsilon) + cosOmegaMoonOnSiniMoon*cos(wpsilon);
    f12=e11*sin(wpsilon) +  ettaAster*cos(wpsilon);
    
    %% Kepler's equation for Moon
    i = 2;
    excentrAnomMoon = qMoon;
    while abs( qMoon+eMoon*sin( excentrAnomMoon ) - excentrAnomMoon ) >= 1e-8
        excentrAnomMoon = qMoon+eMoon*sin( excentrAnomMoon );
        i=i+1;
    end
    
    %% Other constants
    sinFiMoon = ( sqrt(1-eMoon^2)*sin(excentrAnomMoon) )/( 1-eMoon*cos(excentrAnomMoon) );
    cosFiMoon = ( cos(excentrAnomMoon) - eMoon )/( 1-eMoon*cos(excentrAnomMoon) );
    
    eMoon = (sinFiMoon*cos(gMoon)+cosFiMoon*sin(gMoon) )*e11+( cosFiMoon*cos(gMoon) - sinFiMoon*sin(gMoon) )*e12;
    ettaMoon = ( sinFiMoon*cos(gMoon) + cosFiMoon*sin(gMoon) )*etta11+( cosFiMoon*cos(gMoon) - sinFiMoon*sin(gMoon) )*etta12;
    fMoon = ( sinFiMoon*cos(gMoon) + cosFiMoon*sin(gMoon) )*f11+( cosFiMoon*cos(gMoon) - sinFiMoon*sin(gMoon) )*f12;
    rMoon = aMoon*( 1-eMoon*cos(excentrAnomMoon) );
    
    %% Calculation of accelerations
    
    xNormMoon = coordinates( 1 ) / rMoon;
    yNormMoon = coordinates( 2 ) / rMoon;
    zNormMoon = coordinates( 3 ) / rMoon;
    muNormMoon = muMoon / rMoon^2;
    deltaMoon = sqrt( ( eMoon - xNormMoon )^2 + ( ettaMoon - yNormMoon )^2 + ( fMoon - zNormMoon )^2 );
    
    influence(1,1) = muNormMoon*( ( eMoon - xNormMoon )/deltaMoon^3 - eMoon );
    influence(2,1) = muNormMoon*( ( ettaMoon - yNormMoon )/deltaMoon^3 - ettaMoon );
    influence(3,1) = muNormMoon*( ( fMoon - zNormMoon )/deltaMoon^3 - fMoon );
end
