function [ tRel ] = toaRelativeEffects( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory)
    % toaRelativeEffects. Calculate relativistic effects to signal from X-Ray sources.
    %
    %   Allow to consider of relative effect for time of arrival for specific x-ray source. Take into account following effects:
    %       1. are referred to as Roemer delay (the first-order Doppler delay, the effects of annual parallax);
    %       2. are referred to as the Shapiro delay% effect.
    %
    %   [ tRel ] = toaRelativeEffects( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory)
    %
    %   INPUT
    %       xRaySources           array of x-ray sources (every item should be instance of the XRaySource);
    %       earthEphemeris        earth ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       sunEphemeris          sun ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       spaceshipTrajectory   state space vector of spaceship (array of 1-st, 2-nd, 3-d - trajectory coordinate vectors in [km]).
    %
    %   OUTPUT
    %       tRel    array of relative effect for each x-ray source.
    %
    narginchk(4, 4);
    
    capacity  = size(spaceshipTrajectory, 2);
    dimension = length(xRaySources);
    
    tRel = zeros(dimension, capacity);
    for i = 1:dimension
        tRel(i, :) = relativeEffects(xRaySources(i), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
    end
end

function tRel = relativeEffects(x, capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory)
    c         = 299792.458; % speed of light [km/sec]
    earthR    = [earthEphemeris.x earthEphemeris.y earthEphemeris.z]';
    sunR      = [sunEphemeris.x sunEphemeris.y sunEphemeris.z]';
    rSc       = earthR + spaceshipTrajectory;
    
    if capacity == 1
        xNorm = x.Normal';
    else
        xNorm = cvecrep(x.Normal', capacity);
    end
    
    firstPart  = dotProduct(xNorm, rSc, 1).^2 - normOfEveryRow(rSc, 1).^2;
    secondPart = dotProduct(xNorm, sunR, 1).*dotProduct(xNorm, rSc, 1) - dotProduct(sunR, rSc, 1);
    %         firstPart  = dot(xNorm, rSc, 2).^2 - normOfEveryRow(rSc).^2;
    %         secondPart = dot(xNorm, sunR, 2).*dot(xNorm, rSc, 2) - dot(sunR, rSc, 2);
    tRel = 1/(2*c*x.distance)*firstPart + 1/(c*x.distance)*secondPart;
end
