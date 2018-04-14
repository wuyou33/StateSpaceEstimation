function [ normal ] = norm_from_lat_lon(latitude, longitude)
    % norm_from_lat_lon. calculate normal (the line of sight) from longitude and latitude matlab.
    %
    %   [ normal ] = norm_from_lat_lon( latitude, longitude )
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
    [x, y, z] = lat_long_to_cartesian(latitude, longitude);
    normal = [x, y, z];
end
