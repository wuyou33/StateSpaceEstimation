function [ qout ] = quaternionNormalize( q )
%%  QUATNORMALIZE Normalize a quaternion.
%   N = quaternionNormilize( Q ) calculates the normalized quaternion, N, for a
%   given quaternion, Q.  Input Q is an 1-by-4 matrix containing
%   quaternions.  N returns an 1-by-4 matrix of normalized quaternions.
%   Each element of Q must be a real number.  Additionally, Q has its
%   scalar number as the first column.
%%
qMod = norm(q, 2);
qout = q / qMod;