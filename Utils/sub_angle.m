function c = sub_angle(a, b)
    % sub_angle. 
    % Addition function for 'angle space' sigma-points expressed in radians.
    % This needed to deal with the angular discontinuety at +pi and -pi radians.
    %
    c = a - b;
    two_pi = 2*pi;
    
    idx_1 = c > pi;
    idx_2 = c < -pi;
    
    c(idx_1) = c(idx_1) - two_pi;
    c(idx_2) = c(idx_2) + two_pi;
end
