function y = unwrap_angle_from_min_pi_to_pi(x)
    % unwrap_angle_from_min_pi_to_pi. Unwraps the input angle to an angle between -pi and pi degrees.
    %
    %   y = unwrap_angle_from_min_pi_to_pi(x)
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
