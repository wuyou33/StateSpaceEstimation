function [ diffToa ] = calculateDiffToa( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory )
    % calculateDiffToa. Calculate difference between
    %   - time of arrrival to spacecraft
    %   and
    %   - time of arrival to ssb (solar system baricenter)
    %   for every x-ray source.
    %
    %   Calculated by following expression:
    %       dToa = (n*(r + re))/c + tRel,
    %       where 
    %           r       position sc relative to earth,
    %           re      position earth relative to ssb,
    %           n       line of sight to pulsar,
    %           tRel    relative effects,
    %           c       light speed.
    %
    %   [ diffToa ] = calculateDiffToa( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory )
    %
    %   INPUT
    %       xRaySources            array of x-ray sources (every item should be instance of the XRaySource);
    %       earthEphemeris         earth ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       sunEphemeris           sun ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       spaceshipTrajectory    state space vector of spaceship (array of 1-st, 2-nd, 3-d - trajectory coordinate vectors in [km]).
    %
    %   OUTPUT
    %       diffToa     array of differences between toa on spaceship and toa on ssb.
    %
    narginchk(4, 4);
    
    c               = 299792.458; % speed of light [km/sec]
    earthPos        = [earthEphemeris.x, earthEphemeris.y, earthEphemeris.z]';
    dimension       = length(xRaySources);
    capacity        = size(spaceshipTrajectory, 2);
    tRel            = toaRelativeEffects(xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory);
    
    diffToa = zeros(dimension, capacity);
    
    for i = 1:dimension
        x = xRaySources(i);
        xNorm = cvecrep(x.Normal, capacity);
        
        % spaceshipTrajectory - position SC relative to Earth. Hence need to convert to inertial frame (ICRF) with center in SSB.
        % To convert need to: position_vect_sc_to_ssb = position_vect_earth_to_ssb + position_vect_sc_to_earth
        diffToa(i, :) = (dot(xNorm, earthPos + spaceshipTrajectory(1:3, :), 1)) / c + tRel(i, :);
    end
end
