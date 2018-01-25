function [ influence ] = MoonInfluence( time, coordinates )
    % The function calculates the Moon influence to the satellite trajectory
    % Constant
    muMoon = 4902.835; % [km^3/sec^2] - constant of the gravitational field of the Moon
    aMoon = 3.84385243e5; % [] the semimajor axis of the Moon's orbit
    eMoon = 0.054900489; % [] eccentricity of the Moon orbit
    iMoon = 0.0898041080; % [rad] mean inclination of the Moon's orbit to the plane of the ecliptic
    
    qMoon = 2.3555557435 + 8328.6914257190*time + 0.0001545547*time^2; % [rad] the average anomaly of the Moon
    omegaMoon = 2.1824391966 - 33.7570459536*time + 0.0000362262*time^2; % [rad] the average longitude of the ascending node of the Moon
    gMoon = 1.4547885346 + 71.0176852437*time - 0.0001801481*time^2; % [rad] the longitude of the perigee of the Moon's orbit
    
    epsilon = 0.4090926006 - 0.0002270711*time; % [rad] the average inclination of the ecliptic to the equator
    
    cosOmegaMoonOnSiniMoon = cos(omegaMoon) * sin(iMoon);
    ettaAster = sin(omegaMoon) * sin(iMoon);
    eAster = cos(iMoon);
    
    e11 = sin(omegaMoon) * cos(omegaMoon) * ( 1-cos(iMoon) );
    e12 = 1 - ( sin(omegaMoon) )^2 * ( 1-cos(iMoon) );
    
    etta11=eAster*cos(epsilon) - cosOmegaMoonOnSiniMoon * sin(epsilon);
    etta12=e11*cos(epsilon) + ettaAster*sin(epsilon);
    
    f11=eAster*sin(epsilon) + cosOmegaMoonOnSiniMoon * cos(epsilon);
    f12=e11*sin(epsilon) +  ettaAster*cos(epsilon);
    
    %% Kepler's equation for Moon
    i = 2;
    excentrAnomMoon = qMoon;
    while abs( qMoon+eMoon*sin(excentrAnomMoon) - excentrAnomMoon ) >= 1e-8
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
    
    influence(1, 1) = muNormMoon*( ( eMoon - xNormMoon )/deltaMoon^3 - eMoon );
    influence(2, 1) = muNormMoon*( ( ettaMoon - yNormMoon )/deltaMoon^3 - ettaMoon );
    influence(3, 1) = muNormMoon*( ( fMoon - zNormMoon )/deltaMoon^3 - fMoon );
end
