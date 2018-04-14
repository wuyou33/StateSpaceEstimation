function x_ray_cov = x_ray_dtoa_covariance(x_ray_sources,  detector_area, time_bucket, background_photon_rate, error_budget)
    % xRayToaCovariance. Calcualate TOA covariance for each X-Ray source (covariance of every source is independent).
    %   
    %   Reference: Noise Analysis for X-ray Navigation Systems. John Hanson; Suneel Sheikh; Paul Graven.
    %       2008 IEEE/ION Position, Location and Navigation Symposium. 5-8 May 2008. pp. 704-713.
    %
    %   Calculate covariance for measurement of time of arrival for array of some specific x-Ray source (pulsar or quasar).
    %   Covariance (cov) of time of arriavel calculated by following expression:
    %       cov = gs^2 * T / (A*dt*s) + gb^2 * T / (A*dt*s) * b/s
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
    %   x_ray_cov = xRayToaCovariance(x_ray_sources,  detector_area, time_bucket, background_photon_rate, error_budget)
    %
    %   INPUT
    %       x_ray_sources               array of the instances of the <X_RaySource> (every object describe x-Ray source parameter);
    %       detector_area               is the detector area in [m*m];
    %       time_bucket                 t is the length of the observation in [sec];
    %       background_photon_rate      is the background rate (including detector noise, the diffuse X-ray background,
    %                                       un-cancelled cosmic ray events and steady emission from the pulsar, bc) [photons/cm^2/sec];
    %       error_budget                error budget, allow to model errors from sources, which are not considered in main model [%].
    %
    %   OUTPUT
    %       x_ray_cov   array of noise covariance for each x-ray source.
    %
    narginchk(4, 5);
    
    if nargin == 4
        error_budget = 0;
    end
    
    if (time_bucket <= 0)
        error('[ xRayToaNoise::time_bucket ] should be positive integer');
    end
    
    % detectorArea * timeBucket * convertion_meter2_to_cm2 (required coz intensity proportional cm^2)
    a_dt = time_bucket * detector_area * 1e4;
    source_count = length(x_ray_sources);
    
    covariance_row = zeros(source_count, 1);
    for i = 1:source_count
        x_ray_source = x_ray_sources(i);
        source_part = x_ray_source.gSource^2 * x_ray_source.period / (a_dt * x_ray_source.intensity);
        background_part = x_ray_source.gBackgr^2 * x_ray_source.period / (a_dt * x_ray_source.intensity) * background_photon_rate/x_ray_source.intensity;
        covariance_row(i) = (1 + error_budget / 100)^2 * (source_part + background_part);
    end
    
    % Under the assumption that observations of every X-Ray source are independent
    x_ray_cov = diag(covariance_row);
end
