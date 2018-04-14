classdef PhotonDetector < handle
    % Describe photon detector for X-Ray navigation.
    % Source can be Pulsar or Quasar (with high stable EM emitting)
    
    properties(SetAccess = private)
        x_ray_sources;
        detector_area;
        time_bucket;
        background_photon_rate;
        error_budget;
        dtoa_cov;
    end
    
    methods(Access = public)
        function obj = PhotonDetector(args)
            % Create a new instance of the XRayDetector
            %   args:
            %       x_ray_sources               array of x-ray sources;
            %       detector_area               detector area (cm^2);
            %       time_bucket                 size of the bins used to count photons, ie total observed time (sec);
            %       background_photon_rate      is the background rate (including detector noise, the diffuse X-ray background,
            %                                       un-cancelled cosmic ray events and steady emission from the pulsar, bc) (photons/cm^2/sec);
            %       error_budget                error budget.
            %
            
            obj.x_ray_sources           = args.x_ray_sources;
            obj.detector_area           = args.detector_area;
            obj.time_bucket             = args.time_bucket;
            obj.background_photon_rate  = args.background_photon_rate;
            obj.error_budget            = args.error_budget;
            
            obj.dtoa_cov = x_ray_dtoa_covariance(args.x_ray_sources, args.detector_area, args.time_bucket, args.background_photon_rate, args.error_budget);
        end
        
        function dtoa = eval_dtoa(this, trajectory, earth_ephemeris, sun_ephemeris)
            % calculate difference time of arrival (DTOA)
            dtoa = diff_time_of_arrival(this.x_ray_sources, earth_ephemeris, sun_ephemeris, trajectory(1:3, :));
        end
        
        function noise = generate_noise(this, capacity)
            source_count = length(this.x_ray_sources);
            noise = chol(this.dtoa_cov, 'lower') * randn(source_count, capacity);
        end
    end
    
end
