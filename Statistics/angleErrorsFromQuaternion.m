function [ angleErros ] = angleErrorsFromQuaternion( quaternionMeasurement, quaternionTrue )
    % angleErrorsFromQuaternion. Calculate angle (yaw, pitch, roll) errors between estimated quaternion and true quaternion.
    %
    %   [ angleErros ] = angleErrorsFromQuaternion( quaternionMeasurement, quaternionTrue )
    %
    %   INPUT
    %       quaternionMeasurement   measurement quaternion;
    %       quaternionTrue          true quaternion.
    %
    %   OUTPUT
    %      angleErros   array of angle (yaw, pitch, roll) errors.
    %
    [yaw, pitch, roll]          = quat2angle(quaternionTrue');
    [yawEst, pitchEst, rollEst] = quat2angle(quaternionMeasurement');
    
    angleErros = ([yaw, pitch, roll] - [yawEst, pitchEst, rollEst])';
end
