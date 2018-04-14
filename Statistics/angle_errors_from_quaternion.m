function [ angle_erros ] = angle_errors_from_quaternion( quaternion_estimated, quaternion_true )
    % angle_errors_from_quaternion. Calculate angle (yaw, pitch, roll) errors between estimated quaternion and true quaternion.
    %
    %   [ angle_erros ] = angle_errors_from_quaternion( quaternion_estimated, quaternion_true )
    %
    %   INPUT
    %       quaternion_estimated   measurement quaternion;
    %       quaternion_true          true quaternion.
    %
    %   OUTPUT
    %      angle_erros   array of angle (yaw, pitch, roll) errors.
    %
    [yaw, pitch, roll]          = quat2angle(quaternion_true');
    [yawEst, pitchEst, rollEst] = quat2angle(quaternion_estimated');
    
    angle_erros = ([yaw, pitch, roll] - [yawEst, pitchEst, rollEst])';
end
