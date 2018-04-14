function [ iinsProcCov ] = getIntegratedInsProcCov( estimatorType, accBiasSigma, gyroBiasSigma )
    % getIntegratedInsCov. Build covariance of process noise of integrated INS navigation system based on estimator (filter) type.
    %
    %   [ iinsProcCov ] = getIntegratedInsProcCov( estimatorType )
    %
    %   INPUT
    %       estimatorType       time of esimation / filtration algorithm;
    %       accBiasSigma        rms of accelaration bias noise;
    %       gyroBiasSigma       rms of gyro bias noise.
    %
    %   OUTPUT
    %       iinsProcCov     covariance matrix.
    %
    narginchk(3, 3);
    
    if string_match(estimatorType, {'ukf', 'cdkf', 'srcdkf', 'ckf', 'sckf', 'kf', 'srukf', 'cqkf', 'fdckf', 'sghqf'})
        iinsProcCov = [(1e-2*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (7e-4*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (1e-5*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) 1e2*diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) 1e2*diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif string_match(estimatorType, {'ekf'})
        iinsProcCov = [(1e0*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (1e-3*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (1e-5*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) 1e2*diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) 1e2*diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-10*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-10*eye(3)).^2 ... gyro scale factor [-]
            ];
    elseif string_match(estimatorType, {'pf', 'sppf', 'gspf', 'gmsppf'})
        iinsProcCov = [(1.1e-1*eye(3)).^2 zeros(3, 19); ... distance error [km]^2
            zeros(3, 3) (1e-3*eye(3)).^2 zeros(3, 16); ... velocity error [km/sec]^2
            zeros(4, 6) (1e-6*eye(4)).^2 zeros(4, 12); ... quaternion error [-]
            zeros(3, 10) diag(accBiasSigma).^2 zeros(3, 9); ... acceler bias [km/sec^2]^2
            zeros(3, 13) diag(gyroBiasSigma).^2 zeros(3, 6); ... gyro bias [rad/sec]^2
            zeros(3, 16) (1e-15*eye(3)).^2 zeros(3, 3); ... acceler scale factor [-]
            zeros(3, 19) (1e-15*eye(3)).^2 ... gyro scale factor [-]
            ];
    else
        error('[ getIntegratedInsProcCov :: estimatorType ]. Unknown filtration algorithm');
    end
end
