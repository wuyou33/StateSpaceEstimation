function y = ang180(x)
    % ang180. Unwraps the input angle to an angle between -pi and pi degrees.
    %
    %   y = ang180(x)
    %
    %   INPUT
    %       x   angle.
    %
    %   OUTPUT
    %       y   unwrapped angle between -pi and pi.
    %
    y = mod(x-(1e-10)-pi, 2*pi) - pi + (1e-10);
    % The term 1e-10 is to set an input of -180 to 180.
end
