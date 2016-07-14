function [ matrix ] = quaternion2BMatrix( q )
    %QUATERNION2BMATRIX Convert quaternion q to special B matrix (required for INS SNS integrated system).
   matrix = 0.5*[-q(2) -q(3) -q(4); q(1) -q(4) q(3); q(4) q(1) -q(2); -q(3) q(2) q(1)];    
end
