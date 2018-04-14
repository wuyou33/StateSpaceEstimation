function [ dtoa ] = diff_time_of_arrival( x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory )
    % diff_time_of_arrival. Calculate difference (DTOA) between
    %   - time of arrrival to spacecraft
    %   and
    %   - time of arrival to ssb (solar system baricenter)
    %   for every x-ray source.
    %
    %   The DTOA calculated by following expression from 
    %       'THE USE OF VARIABLE CELESTIAL X-RAY SOURCES FOR SPACECRAFT NAVIGATION' formula 4.32 (page 149 (175)) Sheikh. Ph.D:
    %       
    %       dtoa = (n*(r_sc_e + re))/c + tRel,
    %       where 
    %           r_sc_e  position of spacecraft relative to the Earth,
    %           re      position of the Earth relative to SSB,
    %           n       line of sight (unit direction) to pulsar,
    %           tRel    relative effects,
    %           c       light speed.
    %
    %   [ dtoa ] = diff_time_of_arrival( x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory )
    %
    %   INPUT
    %       x_ray_sources       array of x-ray sources (every item should be instance of the X_RaySource);
    %       earth_ephemeris     earth ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       sun_ephemeris       sun ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
    %       trajectory          state space vector of spaceship (array of 1-st, 2-nd, 3-d - trajectory coordinate vectors in [km]).
    %
    %   OUTPUT
    %       diffToa     array of differences between toa on spaceship and toa on ssb.
    %
    narginchk(4, 4);
    
    c               = 299792.458; % speed of light [km/sec]
    earth_pos       = [earth_ephemeris.x, earth_ephemeris.y, earth_ephemeris.z]';
    dimension       = length(x_ray_sources);
    capacity        = size(trajectory, 2);
    tRel            = dtoa_relative_effects(x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory);
    
    dtoa = zeros(dimension, capacity);
    
    for i = 1:dimension
        x = x_ray_sources(i);
        x_norm = column_vector_replicate(x.Normal, capacity);
        
        % trajectory - position SC relative to Earth. Hence need to convert to inertial frame (ICRF) with center in SSB.
        % To convert need to (sum of vectors): position_vect_sc_to_ssb = position_vect_earth_to_ssb + position_vect_sc_to_earth
        dtoa(i, :) = (dot(x_norm, earth_pos + trajectory(1:3, :), 1)) / c + tRel(i, :);
    end
end
