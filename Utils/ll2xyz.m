function [ x, y, z ] = ll2xyz( latitude, longitude )
% LL2XYZ Convert latitude and longitude to xyz (assume unit radius)
% INPUT:
%       latitude  - latitude [rad];
%       longitude - longitude [rad];
% OUTPUT:
%       x - x coordinate
%       y - y coordinate
%       z - z coordinate
%% 
    x = sin(latitude) * cos(longitude);
    y = sin(latitude) * sin(longitude);
    z = cos(latitude);
end