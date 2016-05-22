function [ tRel ] = toaRelativeEffects( xRaySources, earthEphemeris, sunEphemeris, spaceshipTrajectory)
% toaRelativeEffects: allow to consider of relative effect for time of
% arrival for specific x-ray source (
% 1. are referred to as Roemer delay (
% the first-order Doppler delay, the effects of annual parallax), 
% 2.are referred to as the Shapiro delay% effect.)
%   
% INPUT:
%       xRaySources         - array of x-ray sources (every item should be instance of the XRaySource);
%       earthEphemeris      - earth ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
%       sunEphemeris        - sun ephemeris (x, y, z - vectors in [km], vx, vy, vz - vectors in [km/sec]);
%       spaceshipTrajectory - state space vector of spaceship (array of 1-st, 2-nd, 3-d - trajectory coordinate vectors in [km]);
% OUTPUT:
%       tRel                - array of relative effect for each x-ray source;
%%
    if (nargin ~= 4); error('[ toaRelativeEffects ] incorrect number of input arg. Should be 3'); end
%%      
        
    [capacity, ~]  = size(spaceshipTrajectory);
    dimension = length(xRaySources);
    
    if dimension == 7 % optimization
        tRel(:, 1) = relativeEffects(xRaySources(1), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 2) = relativeEffects(xRaySources(2), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 3) = relativeEffects(xRaySources(3), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 4) = relativeEffects(xRaySources(4), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 5) = relativeEffects(xRaySources(5), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 6) = relativeEffects(xRaySources(6), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        tRel(:, 7) = relativeEffects(xRaySources(7), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
    else
        tRel      = zeros(capacity, dimension);
        for i = 1:dimension        
            tRel(:, i) = relativeEffects(xRaySources(i), capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory);
        end
    end
end

function tRel = relativeEffects(x, capacity, earthEphemeris, sunEphemeris, spaceshipTrajectory)
    c         = 299792.458; % speed of light [km/sec]   
    earthR    = [earthEphemeris.x earthEphemeris.y earthEphemeris.z];
    sunR      = [sunEphemeris.x sunEphemeris.y sunEphemeris.z];    
    rSc       = earthR + spaceshipTrajectory;    
    
    if capacity == 1
        xNorm = x.Normal;
    else
        xNorm = repmat(x.Normal, capacity, 1);
    end
    
    firstPart  = dotProduct(xNorm, rSc, 2).^2 - normOfEveryRow(rSc).^2;
    secondPart = dotProduct(xNorm, sunR, 2).*dotProduct(xNorm, rSc, 2) - dotProduct(sunR, rSc, 2);
%         firstPart  = dot(xNorm, rSc, 2).^2 - normOfEveryRow(rSc).^2;
%         secondPart = dot(xNorm, sunR, 2).*dot(xNorm, rSc, 2) - dot(sunR, rSc, 2);
    tRel = 1/(2*c*x.distance)*firstPart + 1/(c*x.distance)*secondPart;
end