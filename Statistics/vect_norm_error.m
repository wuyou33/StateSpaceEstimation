function [ error ] = vect_norm_error(measurement, true_vector, multiplier)
    % vect_norm_error. Calculate norm in N-dimension Euclidean space of difference between array of vectors of estimation and array of vectors with true values.
    %
    %   [ error ] = vect_norm_error(measurement, true_vector, multiplier).
    %
    %   INPUT
    %       measurement 	measurement of vector in N dimensil space, where N - is a number of rows. Every column represent measurement;
    %       true_vector  	true values measurement vector (dimension should be same as measurement);
    %       multiplier      <<optional>> scale of measurement, default value is 1.
    %
    %   OUTPUT
    %       error   vector of error of Eucledian norm in N dimensial space (where N - is a number of rows in input matrices).
    %
    narginchk(2, 3);
    
    if nargin < 3; multiplier = 1; end
    
    error = sum( (multiplier*(true_vector - measurement)).^2, 1 ).^0.5;
end
