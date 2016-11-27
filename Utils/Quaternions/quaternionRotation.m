function [ vect ] = quaternionRotation( q, r, direction )
    %  quaternionRotation. Rotate a vector by a quaternion.
    %
    %   [ vect ] = quaternionRotation( q, r, direction )
    %
    %   Calculates the rotated vector, vect, for a quaternion, q, and a vector, r.  q is either a single 1-by4 quaternion. r is either a single 1-by-3 vector.
    %   vect returns an 1-by-3 matrix of rotated vectors. Each element of q and r must be a real number.
    %
    %   INPUT
    %       q           quaternion;
    %       r           3D vector;
    %       direction   determine that row or column vector used (1 means that row vector used; 2 means that column vector used).
    %
    %   OUTPUT
    %       vect    rotated vector.
    %
    %   Examples:
    %
    %      q = [1 0 1 0];
    %      r = [1 1 1];
    %
    directionCosineMatrix = quaternion2DirectionCosineMatrix(q);
    
    if nargin == 2 || direction == 1
        vect = (directionCosineMatrix * r')';
    elseif nargin == 3 && direction == 2
        vect = directionCosineMatrix * r;
    end
end
