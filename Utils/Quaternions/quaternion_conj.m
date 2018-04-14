function [ conj ] = quaternion_conj( q, direction )
    % quaternion_conj. Calculate the conjugate of a quaternion.
    %
    %   [ conj ] = quaternion_conj( q, direction )
    %
    %   Calculates the conjugate, conj, for a given quaternion, q. Input q is an M-by-4 matrix containing M quaternions.
    %   conj returns an M-by-4 matrix of conjugates.  Each element of Q must be a real number. Additionally, q has its scalar number as the first column.
    %
    %   INPUT
    %       q           quaternion;
    %       direction 	determine that row or column is values of quaternion.
    %
    %   OUTPUT
    %       conj    calculated conjugate matrix.
    %
    if nargin == 1 || direction == 1
        conj = [q(:, 1) -q(:, 2:4)];
    elseif direction == 2
        conj = [q(1, :); -q(2:4, :)];
    else
        error('[ quaternion_conj::direction ] unknown value of direction. Should be 1 or 2 or not specified');
    end
end