close all; clc; clearvars; clear memoize; % clear memoize required for memoization
addpath(genpath('./'));

date.day                = 17;
date.mon                = 11;
date.year               = 2017;
time_start              = '00:00:00.000';
time_end                = '05:00:00.000';
time_data_x_ray         = TimeExt(time_start, time_end, 100, date, 1e7); % change refreshSunMoonInfluenceTime to real number
number_of_runs          = 4;
second_in_one_minute    = 4;
esitimated_params       = 2;
draw_iterations         = 0;
datapoints_count    	= time_data_x_ray.SimulationNumber;
err_arr                 = zeros(number_of_runs, esitimated_params, datapoints_count - 1);
mass                    = 200; % [kg]

filter_types = {'ukf'};
filter_type = filter_types{1};

b_det   = 0.1; % Detector Background Rate. [photon*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^-2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^-2*sec^-1]
x_ray_source_count      = 4;
time_bucket             = 1e5; % [sec]
detector_area           = 1; % [m^2]
background_photon_rate  = b_det / (detector_area*1e4) + b_diff + b_cosm; % detector_area have to be converted to cm*cm
error_budget = 0; % [%]

gravity_model = fit_solar_system_gravity_model(time_data_x_ray);

earth_ephemeris = load_ephemeris('earth', time_data_x_ray, second_in_one_minute/time_data_x_ray.SampleTime);
sun_ephemeris   = load_ephemeris('sun', time_data_x_ray, second_in_one_minute/time_data_x_ray.SampleTime);
x_ray_sources   = load_x_ray_sources(x_ray_source_count);

dtoa_cov = x_ray_dtoa_covariance(x_ray_sources, detector_area, time_bucket, background_photon_rate);

detector_args.x_ray_sources = x_ray_sources;
detector_args.detector_area = detector_area;
detector_args.time_bucket = time_bucket;
detector_args.background_photon_rate = background_photon_rate;
detector_args.error_budget = error_budget;

detector = PhotonDetector(detector_args);

x_ray_init_arg.mass             = mass;
x_ray_init_arg.sample_time      = time_data_x_ray.SampleTime;
x_ray_init_arg.gravity_model    = gravity_model;
x_ray_init_arg.photon_detector  = detector;
x_ray_init_arg.earth_ephemeris  = earth_ephemeris;
x_ray_init_arg.sun_ephemeris    = sun_ephemeris;

x_ray_init_arg.current_epoch = time_data_x_ray.get_current_epoch(time_data_x_ray.Time(1));
x_ray_init_arg.current_time  = time_data_x_ray.Time(1);

x_nav_model = gssm_x_nav('init', x_ray_init_arg);

inference_model_args.type  = 'state';
inference_model_args.tag   = 'State estimation for X-Ray navigation system';
inference_model_args.model = x_nav_model;

inference_model = inference_model_generator(inference_model_args);
[x_noise_infer, z_noise_infer, inference_model] = inference_noise_generator(inference_model, filter_type);

alpha_list = 0:6;
beta_list = 1:6;
kappa_list = [0.001; 0.01; 0.1; 1];

x_cov_est = [(1e0*eye(3)).^2 zeros(3); zeros(3) (1e-5*eye(3)).^2]; % [km^2  n/a; n/a (km/sec)^2];

for alpha_idx = 1:length(alpha_list)
    for beta_idx = 1:length(beta_list)
        for kappa_idx = 1:length(kappa_list)
            filter_params.alpha = alpha_list(alpha_idx); % scale factor (UKF parameter) 1e-3
            filter_params.beta  = beta_list(beta_idx); % 2 is a optimal setting for Gaussian priors (UKF parameter)
            filter_params.kappa = kappa_list(kappa_idx); % 0 is optimal for state dimension = 2 (UKF parameter)
            
            filter_params.x_cov_initial = x_cov_est;
            
            for i = 1 : number_of_runs
                initial_x_ray_orbit = load_initial_orbit() + [1e-2*randn(3, 1); 1e-5*randn(3, 1); 1e-7*randn(4, 1)];
                initial_x_ray_orbit = initial_x_ray_orbit(1:6);
                
                x_nav_imitator = X_RayNavImitator(x_nav_model);
                [x, z] = x_nav_imitator.resolve(initial_x_ray_orbit, time_data_x_ray);
                x_init = x(:, 1) + sqrt(diag(x_cov_est)) .* randn(inference_model.stateDimension, 1);
                
                % Call inference algorithm / estimator
                x_ray_nav_system = X_RayNavigationSystem.create(filter_type, inference_model, x_noise_infer, z_noise_infer, time_data_x_ray);
                x_est = x_ray_nav_system.resolve(x_init, z, filter_params);
                
                % |r| error, ie norm (module) of attitute vector error, [m]
                err_arr(i, 1, :) = 1e3*(sqrt(sum((x_est(1:3, 2 : end)).^2, 1)) - sqrt(sum((x(1:3, 2 : end)).^2, 1)));
                
                % |v| error, ie (module) of velocity vector error, [m / sec]
                err_arr(i, 2, :) = 1e3*(sqrt(sum((x_est(4:6, 2 : end)).^2, 1)) - sqrt(sum((x(4:6, 2 : end)).^2, 1)));
            end
            
            rmse_x_err_attitude = sqrt( mean( ( squeeze(err_arr(:, 1, :)) ).^2) );
            fprintf('alpha = %d; beta = %d; kappa = %d; RMSE attitude = %d\n', filter_params.alpha, filter_params.beta, filter_params.kappa, rmse_x_err_attitude(end));
            
            rmse_x_err_velocity = sqrt( mean( ( squeeze(err_arr(:, 2, :)) ).^2) );
            fprintf('alpha = %d; beta = %d; kappa = %d; RMSE velocity = %d\n', filter_params.alpha, filter_params.beta, filter_params.kappa, rmse_x_err_velocity(end));
        end
    end
end