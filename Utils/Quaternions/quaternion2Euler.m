function [ eul ] = quaternion2Euler( q, sequence )
    % eul Convert quaternion to Euler angles
    %   [ eul ] = quaternion2Euler( q, sequence ). Converts a unit quaternion rotation into the corresponding Euler angles by the axis rotation sequence.
    %   The input, Q, is an 4-by-1 matrix containing 1 quaternion.
    %   Quaternion represents a 3D rotation and is of the form q = [w x y z],
    %   with a scalar number as the first value. Each element of Q must be a real number.
    %   The output, EUL, is an 3-by-1 array of Euler rotation angles with each
    %   row representing one Euler angle set. Rotation angles are in radians.
    %
    %   The default rotation sequence is 'ZYX', where the order of rotation
    %   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
    %
    %   The following rotation sequences, SEQ, are supported: 'ZYX' and 'ZYZ'.
    %
    narginchk(1, 2);
    
    if (nargin() == 1)
        sequence = 'ZYX';
    end
    
    % Normalize the quaternions
    q = quaternionNormalize(q);
    
    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);
    
    % Pre-allocate output
    
    % The parsed sequence will be in all upper-case letters and validated
    switch sequence
        case 'ZYX'
            aSinInput = -2*(qx.*qz-qw.*qy);
            aSinInput(aSinInput > 1) = 1;
            
            eul = [ atan2( 2*(qx.*qy+qw.*qz), qw.^2 + qx.^2 - qy.^2 - qz.^2 );...
                asin( aSinInput ); ...
                atan2( 2*(qy.*qz+qw.*qx), qw.^2 - qx.^2 - qy.^2 + qz.^2 )];
            
        case 'ZYZ'
            R = quat2rotm(q');
            eul = rotm2eul(R, 'ZYZ')';
            
        otherwise
            error('[ quaternion2Euler::sequence ] unknown sequence');
    end
end
