function [ ephemeris ] = loadEphemeris(type, capacity, repeatInterval)
% loadEarthEphemeris: load ephemeris (earth or sun) from
% "horizons_results.csv" or "horizons_results_sun.csv"
%   For detailed explanation about ephemeris see "horizons_results.txt" or
%   horizons_results_sun.csv
%   Ephemeris generated for 1 day with interval - 1 minute;
% INPUT:
%       capacity         - size of required data;
%       repeatInterval   - repeat interval of each vector in ephemeris (used for data extrapolation from minutes to required time interval);
% OUTPUT:
%       1 row - x [km]
%       2 row - y [km]
%       3 row - z [km]
%       4 row - vx [km/s]
%       5 row - vy [km/s]
%       6 row - vz [km/s]
%%
    switch type
        case 'earth'
            path = 'horizons_results.csv';
        case 'sun'
            path = 'horizons_results_sun.csv';
        otherwise
            error('[ loadEphemeris ] not supported type');
    end
     
    data = csvread(path);
    [~, columnCount] = size(data);
    % data generated in http://ssd.jpl.nasa.gov/horizons.cgi in csv has comma after last column, therefore matlab think that it's another new empty column
    columnCount = columnCount - 1; 
    
    extrapolated = data(:, 1:columnCount);
    
    if (nargin == 3 && repeatInterval > 0)
        extrapolated = repelem(extrapolated, repeatInterval, 1); 
    end
    
    ephemeris.x  = extrapolated(1:capacity, 1);
    ephemeris.y  = extrapolated(1:capacity, 2);
    ephemeris.z  = extrapolated(1:capacity, 3);    
    ephemeris.vx = extrapolated(1:capacity, 4);
    ephemeris.vy = extrapolated(1:capacity, 5);
    ephemeris.vz = extrapolated(1:capacity, 6);
end