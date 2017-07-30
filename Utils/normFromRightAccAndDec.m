function [ normal ] = normFromRightAccAndDec( rightAscension, declination )
    % normFromRightAccAndDec. calculate normal (the line of sight) from longitude and latitude matlab.
    %
    %   [ normal ] = normFromLatLon( latitude, longitude )
    %
    %   INPUT
    %       rightAscension  	right ascension angle [rad];
    %       declination         declination angle [rad].
    %
    %   OUTPUT
    %       normal      vector [x, y, z] with following coordinates:
    %           x 	x coordinate;
    %           y 	y coordinate;
    %           z 	z coordinate.
    %
    
    normal = [...
        cos(declination) * cos(rightAscension); ...
        cos(declination) * sin(rightAscension); ...
        sin(declination)];
end
