function [ normal ] = norm_from_right_acc_and_dec(right_ascension, declination)
    % norm_from_right_acc_and_dec. calculate normal (the line of sight) from longitude and latitude matlab.
    %
    %   [ normal ] = norm_from_right_acc_and_dec( latitude, longitude )
    %
    %   INPUT
    %       right_ascension  	right ascension angle [rad];
    %       declination         declination angle [rad].
    %
    %   OUTPUT
    %       normal      vector [x, y, z] with following coordinates:
    %           x 	x coordinate;
    %           y 	y coordinate;
    %           z 	z coordinate.
    %
    
    normal = [...
        cos(declination) * cos(right_ascension); ...
        cos(declination) * sin(right_ascension); ...
        sin(declination)];
end
