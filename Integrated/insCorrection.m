function [ corrected ] = insCorrection( insMeasurement, correctionParams)
%   Correction measurement of strapdown inertial navigation sistem (INS)
%   according with correction parameters
%   
%   INPUT
%       insMeasurement:     measurement of INS
%                               insMeasurement(1:3) - position
%                               insMeasurement(4:6) - velocity
%                               insMeasurement(7:10) - quaternion
%       correctionParams:   correction parameters
%                               correctionParams(1:3) - correction param for position
%                               correctionParams(4:6) - correction param for velocity
%                               correctionParams(7:10) - correction param for quaternion
%   OUTPUT
%       corrected:          corrected measurement of INS
%                               corrected(1:3) - position
%                               corrected(4:6) - velocity
%                               corrected(7:10) - quaternion
%%
    corrected(1:3)  = insMeasurement(1:3) - correctionParams(1:3);
    corrected(4:6)  = insMeasurement(4:6) - correctionParams(4:6);
    corrected(7:10) = quaternionNormalize(quaternionMultiply(insMeasurement(7:10), correctionParams(7:10)));    
end

