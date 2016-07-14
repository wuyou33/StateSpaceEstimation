function [ conj ] = quaternionConj( q, direction )
% quaternionConj: Calculate the conjugate of a quaternion.
%   N = quaternionConj( Q ) calculates the conjugate, N, for a given quaternion, Q.  
%   Input Q is an M-by-4 matrix containing M quaternions.  N returns an 
%   M-by-4 matrix of conjugates.  Each element of Q must be a real number.  
%   Additionally, Q has its scalar number as the first column.
%
% direction - determine that row or column is values of quaternion.
    
    if nargin == 1 || direction == 1
        conj = [q(:, 1) -q(:, 2:4)];
    elseif direction == 2
        conj = [q(1, :); -q(2:4, :)];
    else
        error(' [ quaternionConj ]: unknown value of direction. Should be 1 or 2 or not specified');
    end
end