function [ covariance ] = xRayToaCovariance( xRaySources,  detectorArea, timeBucket, backgroundPhotonRate, errorBudget )
    % xRayToaCovariance. Calcualate TOA covariance for each X-Ray source (covariance of every source is independent).
    %
    %   Calculate covariance for measurement of time of arrival for array of some specific x-Ray source (pulsar or quasar).
    %   RMS of time of arriavel calculated by following expression:
    %       rms = sqrt(gs^2 * T / (A*dt*s) + gb^2 * T / (A*dt*s) * b/s)
    %         where
    %         gs - geometric factor which dependent from source (sec^0.5);
    %         A  - detector area (cm^2);
    %         T  - period of signal x-ray source (sec);
    %         dt - size of the bins used to count photons, ie total observed time (sec);
    %         s  - average x-ray source Flux (photons/cm^2/sec);
    %         b  - is the background rate (including detector noise, the diffuse X-ray background,
    %               un-cancelled cosmic ray events and steady emission from the pulsar, bc) (photons/cm^2/sec);
    %         gb - geometric factor which dependent from background emittion (sec^0.5);
    %
    %   [ covariance ] = xRayToaCovariance( xRaySources,  detectorArea, timeBucket, backgroundPhotonRate )
    %
    %   INPUT
    %       xRaySources            array of the instances of the 'XRaySource' (every object describe x-Ray source parameter);
    %       detectorArea           is the detector area in cm*cm;
    %       timeBucket             t is the length of the observation in sec;
    %       capacity               number of sampling which sould be generated;
    %       backgroundPhotonRate    is the background rate (including detector noise, the diffuse X-ray background,
    %                              un-cancelled cosmic ray events and steady emission from the pulsar, bc) [photons/cm^2/sec];
    %       errorBudget            error budget, allow to model errors from sources, which are not considered in main model [%].
    %   OUTPUT
    %       covariance  array of noise covariance for each x-ray source.
    %
    narginchk(4, 5);
    
    if nargin == 4
        errorBudget = 0;
    end
    
    if (timeBucket <= 0)
        error('[ xRayToaNoise::timeBucket ] timeBucket should be positive integer');
    end
    
    % detectorArea * timeBucket * convertionBetweenMeter2AndCm2 (required coz intensity proportional cm^2)
    a_dt = detectorArea*timeBucket*1e4;
    sourceCount = length(xRaySources);
    
    covarianceRow = zeros(sourceCount, 1);
    for i = 1:sourceCount
        xRaySource = xRaySources(i);
        sourcePart = xRaySource.gSource^2 * xRaySource.period / (a_dt * xRaySource.intensity);
        backgroundPart = xRaySource.gBackgr^2 * xRaySource.period / (a_dt * xRaySource.intensity) * backgroundPhotonRate/xRaySource.intensity;
        covarianceRow(i) = (1 + errorBudget / 100)^2 * (sourcePart + backgroundPart);
    end
    
    % Under the assumption that observations of every X-Ray source are independent
    covariance = diag(covarianceRow);
end
