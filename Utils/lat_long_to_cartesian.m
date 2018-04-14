function [ x, y, z ] = lat_long_to_cartesian(latitude, longitude)
    % lat_long_to_cartesian. Convert latitude and longitude to xyz (assume unit radius).
    %
    %   [ x, y, z ] = lat_long_to_cartesian( latitude, longitude )
    %
    %   INPUT
    %       latitude  	latitude [rad];
    %       longitude 	longitude [rad].
    %
    %   OUTPUT
    %       x 	x coordinate;
    %       y 	y coordinate;
    %       z 	z coordinate.
    %
    x = sin(latitude) * cos(longitude);
    y = sin(latitude) * sin(longitude);
    z = cos(latitude);
end
