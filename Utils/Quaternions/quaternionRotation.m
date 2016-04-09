function [ vect ] = quaternionRotation( q, r )
%%  quaternionRotation: Rotate a vector by a quaternion.
%   vect = QUATROTATE( q, r ) calculates the rotated vector, vect, for a
%   quaternion, q, and a vector, r.  Q is either a single 1-by4 quaternion. 
%   r is either a single 1-by-3 vector.  vect returns an
%   1-by-3 matrix of rotated vectors.  Each element of q and r must be a
%   real number.  
%   
%   Examples:
%
%      q = [1 0 1 0];
%      r = [1 1 1];
%

%%
directionCosineMatrix = quaternion2DirectionCosineMatrix(q);
vect = (directionCosineMatrix*r')';

