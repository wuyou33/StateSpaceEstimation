function [ error ] = vectNormError(measurement, trueVector, multiplier)
    % vectNormError. Calculate norm in N-dimension Euclidean space of difference between array of vectors of estimation and array of vectors with true values.
    %
    %   [ error ] = vectNormError(measurement, trueVector, multiplier).
    %
    %   INPUT
    %       measurement 	measurement of vector in N dimensil space, where N - is a number of rows. Every column represent measurement;
    %       trueVector  	true values measurement vector (dimension should be same as measurement);
    %       multiplier      <<optional>> scale of measurement, default value is 1.
    %
    %   OUTPUT
    %       error   vector of error of Eucledian norm in N dimensial space (where N - is a number of rows in input matrices).
    %
    narginchk(2, 3);
    
    if nargin < 3; multiplier = 1; end
    
    error = sum( (multiplier*(trueVector - measurement)).^2, 1 ).^0.5;
end
