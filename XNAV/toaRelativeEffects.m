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
    c         = 299792.458; % speed of light [km/sec]  
    earthR    = [earthEphemeris.x earthEphemeris.y earthEphemeris.z];
    sunR      = [sunEphemeris.x sunEphemeris.y sunEphemeris.z];    
    rSc       = earthR + spaceshipTrajectory;
    dimension = length(xRaySources);
    [capacity, ~]  = size(spaceshipTrajectory);
    
    tRel      = zeros(capacity, dimension);
    
    for i = 1:dimension        
        x = xRaySources(i);
        xNorm = repmat(x.Normal, capacity, 1);        
        
        firstPart  = dot(xNorm, rSc, 2).^2 - normOfEveryRow(rSc).^2;
        secondPart = dot(xNorm, sunR, 2).*dot(xNorm, rSc, 2) - dot(sunR, rSc, 2);
        
        tRel(:, i) = 1/(2*c*x.distance)*firstPart + 1/(c*x.distance)*secondPart;
    end
end