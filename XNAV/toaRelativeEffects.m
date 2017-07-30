function [ tRel ] = toaRelativeEffects( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory)
    % toaRelativeEffects. Calculate relativistic effects to signal from X-Ray sources.
    %
    %   Allow to consider of relative effect for time of arrival for specific x-ray source. Take into account following effects:
    %       1. are referred to as Roemer delay (the first-order Doppler delay, the effects of annual parallax);
    %       2. are referred to as the Shapiro delay effect.
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
    %   Function to calculate relativistic effects.
    %   'THE USE OF VARIABLE CELESTIAL X-RAY SOURCES FOR SPACECRAFT NAVIGATION' formula 4.32. Sheikh. Ph.D
    %   The first term on the right-hand side of equation is the first order Doppler delay, 
    %   and represents the simple geometric time delay between these two locations. 
    %   The second term is due to the effects of parallax. Together these two terms are referred to as Roemer delay. 
    %   The last term is the Sun?s Shapiro delay effect, which is the additional time delay from the curved light ray path due to the Sun?s gravity field
    %   
    c           = 299792.458; % speed of light [km/sec]
    mu_sun      = 1.327124400189E11; % [km^3*s^?2]
    earthR      = [earthEphemeris.x earthEphemeris.y earthEphemeris.z]'; % [km]
    b           = [sunEphemeris.x sunEphemeris.y sunEphemeris.z]'; % [km]
    r           = earthR + spaceshipTrajectory(1:3, :); % [km]
    d0          = x.distance; % [km]
    
    if capacity == 1
        n = x.Normal;
    else
        n = cvecrep(x.Normal, capacity);
    end
    
    tRel = ( dot(n, r, 1).^2 - norm(r, 1).^2 + 2*dot(n, b, 1).*dot(n, r, 1) - 2*dot(b, r, 1) ) / ( 2*c*d0 ) + ...
        2*mu_sun/c^3*log((dot(n, r, 1) + norm(r, 1)) ./ (dot(n, b, 1)+norm(b, 1)) + ones(1, capacity));
end
