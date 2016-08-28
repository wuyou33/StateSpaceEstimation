function [ radialSpeed ] = dopplerShift( starshipState, earthEphemeris, sunEphemeris )
    %DOPPLERSHIFT Calculate radial velocity of starsip relative to Sun, which related with Doppler shift.
    %   In general, the light of the Sun can be imaged with a spectrometer.
    %   Owing to the relative motion between a light source and a moving object, the spectral lines of light are shifted from the original
    %   positions corresponding to their wavelength, which is called the Doppler shift.
    %   Therefore a radial velocity can be measured from a Doppler shift due to the relative motion between the Sun and
    %   a spacecraft with a resonance-scattering spectrometer, whose accuracy is less than 0.01m/s.
    %   The radial speed, |Vsun|, can be modelled by the following equation
    %
    %                 (R + Re + Rs)*(V + Ve + Vs)
    %   |Vsun| = ---------------------------------------,
    %                        |R + Re + Rs|
    %   where:
    %   Vsun - radial speed;
    %   R  - the position vector of the spacecraft with respect to Earth;
    %   Re - the known position vector of the Earth with respect to the SSB;
    %   Rs - the position vectors of the SSB relative to the Sun;
    %   V  - the velocity vector of the spacecraft with respect to Earth;
    %   Ve - the known velocity vector of the Earth with respect to the SSB;
    %   Vs - the velocity vectors of the SSB relative to the Sun.
    %   SSB - solar system baricenter.
    %
    %   INPUT
    %       starshipState  - vector or array of vector dimension 6-by-N position and velocity of starship for every sample from 1 to N [km] and [km / sec];
    %       earthEphemeris - earth ephemeris, vector or array of vector of position and velocity of Earth for every sample from 1 to N [km] and [km / sec];
    %       sunEphemeris   - sun ephemeris, vector or array of vector of position and velocity of Sun for every sample from 1 to N [km] and [km / sec].
    %
    %   OUTPUT
    %       radialSpeed - radial velocity [km / sec];
    %
    narginchk(3, 3);
    
    re = [earthEphemeris.x earthEphemeris.y earthEphemeris.z]';
    ve = [earthEphemeris.vx earthEphemeris.vy earthEphemeris.vz]';
    
    rs = [sunEphemeris.x sunEphemeris.y sunEphemeris.z]';
    vs = [sunEphemeris.vx sunEphemeris.vy sunEphemeris.vz]';
    
    r = starshipState(1:3, :);
    v = starshipState(4:6, :);
    
    rcum = r + re + rs;
    vcum = v + ve + vs;
    
    radialSpeed = dotProduct(rcum, vcum) ./ normOfEveryRow(rcum);
end
