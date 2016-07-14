function [ error ] = rmsd(measurement, trueVector, multiplier)
    % rmsd
    %   [ error ] = rmsd(measurement)
    %   INPUT:
    %       measurement - measurement of vector in N dimensil space, where N - is a number of rows. Every column represent measurement;
    %       trueVector  - true values measurement vector (dimension should be same as measurement).
    %       multiplier  - <optional> scale of measurement, default value is 1;
    %   OUTPUT:
    %       error  - vector of error of Eucledian norm in N dimensial space (where N - is a number of rows in input matrices).
    
    if nargin < 3; multiplier = 1; end
    
    error = sum( (multiplier*(trueVector - measurement)).^2, 1 ).^0.5;    
end
