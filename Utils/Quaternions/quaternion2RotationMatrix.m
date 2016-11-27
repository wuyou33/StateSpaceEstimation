function [ rotationMatrix ] = quaternion2RotationMatrix( q )
    % quaternion2RotationMatrix. Convert quaternion to rotation matrix.
    %
    %   [ rotationMatrix ] = quaternion2RotationMatrix( q )
    %
    %   Converts a unit quaternion, q, into an orthonormal rotation matrix, rotationMatrix. The input, q, is an 1-by-4 matrix containing a
    %   quaternion (q = [w x y z]). Each element of q must be a real number. The output, rotationMatrix, is an 3-by-3 matrix containing rotation matrix.
    %
    %   INPUT
    %       q   quaternion.
    %
    %   OUTPUT
    %       rotationMatrix  orthonormal rotation matrix.
    %
    % Normalize and transpose the quaternions
    q = quaternionNormalize(q)';
    
    % Reshape the quaternions in the depth dimension
    tempR = [...
        1 - 2*(q(3).^2 + q(4).^2),    2*(q(2).*q(3) - q(1).*q(4)),   2*(q(2).*q(4) + q(1).*q(3)), ...
        2*(q(2).*q(3) + q(1).*q(4)),  1 - 2*(q(2).^2 + q(4).^2),     2*(q(3).*q(4) - q(1).*q(2)), ...
        2*(q(2).*q(4) - q(1).*q(3)),  2*(q(3).*q(4) + q(1).*q(2)),   1 - 2*(q(2).^2 + q(3).^2) ...
        ];
    
    rotationMatrix = permute(reshape(tempR, [3, 3, 1]), [2 1 3]);
end

