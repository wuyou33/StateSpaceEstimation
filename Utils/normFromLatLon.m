function [ normal ] = normFromLatLon( latitude, longitude )
    % normFromLatLon. calculate normal (the line of sight) from longitude and latitude matlab.
    %
    %   [ normal ] = normFromLatLon( latitude, longitude )
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
    [x, y, z] = ll2xyz(latitude, longitude);
    normal = [x, y, z];
end
