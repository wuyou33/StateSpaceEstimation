function [ corrected ] = insCorrection( insMeasurement, correctionParams)
    % insCorrection. Correction measurement of strapdown inertial navigation sistem (INS) according with correction parameters (estimated errors of INS).
    %
    %   [ corrected ] = insCorrection( insMeasurement, correctionParams)
    %
    %   INPUT
    %       insMeasurement      measurement of INS
    %                               insMeasurement(1:3) - position
    %                               insMeasurement(4:6) - velocity
    %                               insMeasurement(7:10) - quaternion
    %       correctionParams    correction parameters
    %                               correctionParams(1:3) - correction param for position
    %                               correctionParams(4:6) - correction param for velocity
    %                               correctionParams(7:10) - correction param for quaternion
    %
    %   OUTPUT
    %       corrected:          corrected measurement of INS
    %                               corrected(1:3) - position
    %                               corrected(4:6) - velocity
    %                               corrected(7:10) - quaternion
    %
    corrected(1:3, 1)  = insMeasurement(1:3) - correctionParams(1:3);
    corrected(4:6, 1)  = insMeasurement(4:6) - correctionParams(4:6);
    
    % Quaternion corrections
    % Crassidis. Eq. 7.34 and A.174a. WITH Brute-force normalization
    insQuat = insMeasurement(7:10);
    antm = [0 insQuat(3) -insQuat(2); -insQuat(3) 0 insQuat(1); insQuat(2) -insQuat(1) 0];
    quat = insQuat + 0.5 .* [insQuat(4)*eye(3) + antm; -1.*[insQuat(1) insQuat(2) insQuat(3)]] * quaternion_to_euler(correctionParams(7:10));
    corrected(7:10, 1) = quaternion_normalize(quat);
    corrected(7:10, 1) = quat;
    
    % corrected(7:10, 1) =  quaternion_normalize(quaternion_multiply(insMeasurement(7:10), quatconj(correctionParams(7:10)')'));
    % corrected(7:10, 1) =  quaternion_multiply(correctionParams(7:10), quatconj(insMeasurement(7:10)')');
    % corrected(7:10, 1) = insMeasurement(7:10);
end
