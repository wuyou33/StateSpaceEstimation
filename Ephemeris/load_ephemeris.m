function [ ephemeris ] = load_ephemeris(type, time_data, repeat_interval)
    % loadEarthEphemeris. Load ephemeris (earth or sun) from "horizons_results.csv" or "horizons_results_sun.csv".
    %
    %   For detailed explanation about ephemeris see "horizons_results.txt" or horizons_results_sun.csv
    %
    %   Ephemeris generated for 1 day with interval - 1 minute.
    %
    %   INPUT
    %       time_data          information about simulation time (start time, end time, simulation number, sample time, etc);
    %       repeatInterval     repeat interval of each vector in ephemeris (used for data extrapolation from minutes to required time interval).
    %
    %   OUTPUT
    %       map:
    %           key - time,
    %           value - struct with following fields:
    %               x - attitude by axis x [km];
    %               y - attitude by axis y [km];
    %               z - attitude by axis z [km];
    %       	    vx - velocity by axis x [km/s];
    %               vy - velocity by axis y [km/s];
    %               vz - velocity by axis z [km/s].
    %
    switch type
        case 'earth'
            path = 'horizons_results.csv';
        case 'sun'
            path = 'horizons_results_sun.csv';
        otherwise
            error('[ loadEphemeris::type ] not supported type');
    end
    
    data = csvread(path);
    [~, columns_count] = size(data);
    % data generated in http://ssd.jpl.nasa.gov/horizons.cgi in csv has comma after last column, therefore matlab think that it's another new empty column
    columns_count = columns_count - 1;
    
    capacity = fix(time_data.SimulationNumber);
    
    data = data(:, 1:columns_count);
    
    if (nargin == 3 && repeat_interval > 0)
        n = length(data)*repeat_interval;
        indexes = linspace(1, n, length(data));
        interpolated = zeros(capacity, 6);
        
        for i = 1:6
            interpolated(:, i) = interp1(indexes, data(:, i), 1:capacity, 'spline');
        end
    else
        interpolated = data;
    end
    
    ephemeris = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    for i = 1:capacity
        ephemeris_item = struct('x', interpolated(i, 1), 'y', interpolated(i, 2), 'z', interpolated(i, 3),...
            'vx', interpolated(i, 4), 'vy', interpolated(i, 5), 'vz', interpolated(i, 6));
        ephemeris(time_data.Time(i)) = ephemeris_item;
    end
end
