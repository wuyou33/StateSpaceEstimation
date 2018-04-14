function [ dtoa_rel ] = dtoa_relative_effects(x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory)
    % toaRelativeEffects. Calculate relativistic effects to signal from X-Ray sources.
    %
    %   Allow to consider of relative effect for time of arrival for specific x-ray source. Take into account following effects:
    %       1. are referred to as Roemer delay (the first-order Doppler delay, the effects of annual parallax);
    %       2. are referred to as the Shapiro delay effect.
    %
    %   [ dtoa_rel ] = dtoa_relative_effects( x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory)
    %
    %   INPUT
    %       x_ray_sources       array of x-ray sources (every item should be instance of the <X_RaySource>);
    %       earth_ephemeris     earth ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       sun_ephemeris       sun ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       trajectory          state space vector of spaceship (array of 1-st, 2-nd, 3-d - trajectory coordinate vectors in [km]).
    %
    %   OUTPUT
    %       dtoa_rel    array of relative effect for dtoa (difference time of arrival between spacecraft and SSB) for each x-ray source.
    %
    narginchk(4, 4);
    
    capacity  = size(trajectory, 2);
    dimension = length(x_ray_sources);
    
    dtoa_rel = zeros(dimension, capacity);
    for i = 1:dimension
        dtoa_rel(i, :) = relative_effects(x_ray_sources(i), capacity, earth_ephemeris, sun_ephemeris, trajectory);
    end
end

function dtoa_rel = relative_effects(x_ray, capacity, earth_ephemeris, sun_ephemeris, trajectory)
    %   Function to calculate relativistic effects.
    %
    %   dtoa_rel = relative_effects(x_ray, capacity, earth_ephemeris, sun_ephemeris, trajectory)
    %
    %   'THE USE OF VARIABLE CELESTIAL X-RAY SOURCES FOR SPACECRAFT NAVIGATION' formula 4.32 (page 149 (175)) Sheikh. Ph.D
    %   The first term on the right-hand side of equation is the first order Doppler delay, 
    %       and represents the simple geometric time delay between these two locations. 
    %   The second term is due to the effects of parallax. Together these two terms are referred to as Roemer delay. 
    %   The last term is the Sun's Shapiro delay effect, which is the additional time delay from the curved light ray path due to the Sun's gravity field.
    %   
    c       = 299792.458; % speed of light [km/sec]
    mu_sun  = 1.327124400189E11; % [km^3*s^2]
    r_e     = [earth_ephemeris.x earth_ephemeris.y earth_ephemeris.z]'; % [km]
    b       = [sun_ephemeris.x sun_ephemeris.y sun_ephemeris.z]'; % [km]
    r       = r_e + trajectory(1:3, :); % [km]
    d0      = x_ray.distance; % [km]
    n = column_vector_replicate(x_ray.Normal, capacity); % unit direction to Pulsar [--]
    
    dtoa_rel = ( dot(n, r, 1).^2 - norm(r, 1).^2 + 2*dot(n, b, 1).*dot(n, r, 1) - 2*dot(b, r, 1) ) / ( 2*c*d0 ) + ...
        2*mu_sun/c^3*log( (dot(n, r, 1) + norm(r, 1)) ./ (dot(n, b, 1)+norm(b, 1)) + ones(1, capacity) );
end
