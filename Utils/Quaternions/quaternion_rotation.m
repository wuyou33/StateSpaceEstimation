function [ vect ] = quaternion_rotation( q, r, direction )
    %  quaternion_rotation. Rotate a vector by a quaternion.
    %
    %   [ vect ] = quaternion_rotation( q, r, direction )
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
    dcm = quaternion_to_direction_cosine_matrix(q);
    
    if nargin == 2 || direction == 1
        vect = (dcm * r')';
    elseif nargin == 3 && direction == 2
        vect = dcm * r;
    end
end
