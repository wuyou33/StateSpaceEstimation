function [ matrix ] = quaternion2BMatrix( q )
    % quaternion2BMatrix. Convert quaternion to special matrix representation.
    %
    %   Convert quaternion q to special B matrix (required for INS SNS integrated system). Allow to change quaternion multiplication of 3D vector and quaternion
    %   on matrix multiplication of 3D vector and matrix.
    %
    %   INPUT
    %       q   quaternion.
    %
    %   OUTPUT
    %       matrix      matrix representation of quaternion.
    %
    narginchk(1, 1);
    matrix = 0.5*[-q(2) -q(3) -q(4); q(1) -q(4) q(3); q(4) q(1) -q(2); -q(3) q(2) q(1)];
end
