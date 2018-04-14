close all; clc; clearvars; clear memoize; % clear memoize required for memoization
addpath(genpath('./'));

date.day                = 17;
date.mon                = 11;
date.year               = 2017;
time_start              = '00:00:00.000';
time_end                = '25:00:00.000';
time_data_x_ray         = TimeExt(time_start, time_end, 100, date, 1e7); % change refreshSunMoonInfluenceTime to real number
number_of_runs          = 4;
second_in_one_minute    = 60;
esitimated_params       = 2;
draw_iterations         = 1;
datapoints_count    	= time_data_x_ray.SimulationNumber;
err_arr                 = zeros(number_of_runs, esitimated_params, datapoints_count - 1);
mass                    = 200; % [kg]

%{'ukf', 'srukf', 'cdkf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'ekf', 'pf', 'sppf', 'gspf', 'gmsppf'};
filter_types = {'ukf'};
filter_type = filter_types{1};

b_det   = 0.1; % Detector Background Rate. [photon*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^-2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^-2*sec^-1]
x_ray_source_count      = 4;
time_bucket             = 1e7; % [sec]
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

model = gssm_x_nav('init', x_ray_init_arg);

inference_model_args.type  = 'state';
inference_model_args.tag   = 'State estimation for X-Ray navigation system';
inference_model_args.model = model;

x_nav_inference_model = inference_model_generator(inference_model_args);
[x_noise_infer, z_noise_infer, x_nav_inference_model] = inference_noise_generator(x_nav_inference_model, filter_type);

tic;
for i = 1 : number_of_runs
    % random permutation
    rng(sum(100*clock), 'v5normal');
    rng(sum(100*clock), 'v5uniform');
    
    initial_x_ray_orbit = load_initial_orbit() + [1e-1*randn(3, 1); 1e-4*randn(3, 1); 1e-7*randn(4, 1)];
    initial_x_ray_orbit = initial_x_ray_orbit(1:6);
    
    x_nav_imitator = X_RayNavImitator(model);
    [x, z] = x_nav_imitator.resolve(initial_x_ray_orbit, time_data_x_ray);
    
    x_cov_est = [(5e0*eye(3)).^2 zeros(3); zeros(3) (1e-5*eye(3)).^2]; % [km^2  n/a; n/a (km/sec)^2];
    x_init = x(:, 1) + sqrt(diag(x_cov_est)) .* randn(x_nav_inference_model.stateDimension, 1);
    
    % build filter parameters
    switch filter_type
        case 'ukf'
            filter_params.alpha = 1e-1; % scale factor (UKF parameter) 1e-3  5e-1
            filter_params.beta  = 3; % 2 is a optimal setting for Gaussian priors (UKF parameter)
            filter_params.kappa = 0.5; % 0 is optimal for state dimension = 2 (UKF parameter)
            filter_params.x_cov_initial = x_cov_est;  
        case 'cdkf'
            x_nav_inference_model.spkfParams = sqrt(3);
            
            for j = 2 : datapoints_count
                current_time = time_data_x_ray.Time(j);
                current_epoch = time_data_x_ray.get_current_epoch(current_time);
                
                x_nav_inference_model.model = x_nav_inference_model.model.set_params(x_nav_inference_model.model, current_epoch, current_time);
                
                [x_est(:, j), x_cov_est] = cdkf(x_est(:, j-1), x_cov_est, x_noise_infer, z_noise_infer, z(:, j), x_nav_inference_model, [], []);
            end
        case 'pf'
            particlesCount = 5e4; % 200
            particle_set.particlesNum = particlesCount;
            particle_set.particles = randn(x_nav_inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
            particle_set.weights = column_vector_replicate(1 / particlesCount, particlesCount);
            x_nav_inference_model.resampleThreshold = 0.1;
            x_nav_inference_model.estimateType = 'mean';
            x_nav_inference_model.resampleMethod = 'residual-kitagawa';
            
            for j = 2 : datapoints_count
                [x_est(:, j), particle_set] = pf(particle_set, x_noise_infer, z_noise_infer, z(j), u1(j-1), u2(j), x_nav_inference_model);
            end
        case 'gspf'
            particlesCount = 10000; % 200
            particle_set.particlesNum = particlesCount;
            initial_particles = randn(x_nav_inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
            
            % fit a n component GMM to initial state distribution
            particle_set.stateGMM = gmm_fit(initial_particles, 4, [0.001 10], 'sqrt');
            
            x_nav_inference_model.estimateType = 'mean';
            x_nav_inference_model.threshold = 0.001;
            x_nav_inference_model.resampleMethod = 'residual';
            
            gspf_arg.type = 'gmm';
            gspf_arg.covarianceType = 'sqrt';
            gspf_arg.dimension = model.processNoiseDimension;
            gspf_arg.mixtureCount = 4;
            gspf_arg.mean = column_vector_replicate(model.processNoise.mean, gspf_arg.mixtureCount);
            gspf_arg.covariance = zeros(gspf_arg.dimension, gspf_arg.dimension, gspf_arg.mixtureCount);
            gspf_arg.covariance(:, :, 1) = 2 * model.processNoise.covariance(:, :, 1);
            gspf_arg.covariance(:, :, 2) = 0.5 * model.processNoise.covariance(:, :, 1);
            gspf_arg.covariance(:, :, 3) = 1.5 * model.processNoise.covariance(:, :, 1);
            gspf_arg.covariance(:, :, 4) = 1.25 * model.processNoise.covariance(:, :, 1);
            gspf_arg.weights = [0.45 0.45 0.05 0.05];
            gspf_process_noise = generate_noise_model(gspf_arg);
            
            for j = 2 : datapoints_count
                [x_est(:, j), particle_set] = gspf(particle_set, gspf_process_noise, z_noise_infer, z(j), u1(j-1), u2(j), x_nav_inference_model);
            end
        case 'sppf'
            particlesCount = 400; % 200
            particle_set.particlesNum = particlesCount;
            particle_set.particles    = column_vector_replicate(x_est(:, 1), particlesCount);
            particle_set.particlesCov = repmat(eye(x_nav_inference_model.stateDimension), [1 1 particlesCount]);
            
            proc_noise_gauss.covariance = sqrt(2*3/4);
            observ_noise_gauss.covariance = sqrt(1e-1);
            
            [proc_noise_gauss, observ_noise_gauss, ~] = inference_noise_generator(x_nav_inference_model, 'srukf');
            
            particle_set.processNoise = proc_noise_gauss;
            particle_set.observationNoise = observ_noise_gauss;
            particle_set.weights = column_vector_replicate(1 / particlesCount, particlesCount);
            
            x_nav_inference_model.spkfType = 'srukf';
            x_nav_inference_model.spkfParams = [1 2 0]; % [1 0 2];
            x_nav_inference_model.resampleThreshold = 1;
            x_nav_inference_model.estimateType = 'mean';
            x_nav_inference_model.resampleMethod = 'residual';
            
            for j = 2 : datapoints_count
                [x_est(:, j), particle_set] = sppf(particle_set, x_noise_infer, z_noise_infer, z(j), u1(j-1), u2(j), x_nav_inference_model);
            end
        case 'gmsppf'
            particlesCount = 10000; % 200
            particle_set.particlesNum = particlesCount;
            
            initial_particles = randn(x_nav_inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
            
            tempCov = zeros(1, 1, 2);
            tempCov(:, :, 1) = sqrt(2);
            tempCov(:, :, 2) = 1;
            
            particle_set.stateGMM = gmm_fit(initial_particles, 4, [0.001 10], 'sqrt');
            
            x_nav_inference_model.estimateType = 'mean';
            x_nav_inference_model.spkfType = 'srukf';
            x_nav_inference_model.spkfParams = [1 2 0]; % [1 0 2];
            x_nav_inference_model.resampleMethod = 'residual';
            
            gmsppf_arg.type = 'gmm';
            gmsppf_arg.covarianceType = 'sqrt';
            gmsppf_arg.dimension = model.processNoiseDimension;
            gmsppf_arg.mixtureCount = 4;
            gmsppf_arg.mean = column_vector_replicate(model.processNoise.mean, gmsppf_arg.mixtureCount);
            gmsppf_arg.covariance = zeros(gmsppf_arg.dimension, gmsppf_arg.dimension, gmsppf_arg.mixtureCount);
            gmsppf_arg.covariance(:, :, 1) = 2 * model.processNoise.covariance(:, :, 1);
            gmsppf_arg.covariance(:, :, 2) = 0.5 * model.processNoise.covariance(:, :, 1);
            gmsppf_arg.covariance(:, :, 3) = 1.5 * model.processNoise.covariance(:, :, 1);
            gmsppf_arg.covariance(:, :, 4) = 1.25 * model.processNoise.covariance(:, :, 1);
            gmsppf_arg.weights = [0.45 0.45 0.05 0.05];
            gmsppf_process_noise = generate_noise_model(gmsppf_arg);
            
            for j = 2 : datapoints_count
                [x_est(:, j), particle_set] = gmsppf(particle_set, gmsppf_process_noise, z_noise_infer, z(j), u1(j-1), u2(j), x_nav_inference_model);
            end
        otherwise
            error('[ filter_type ] Unknown / not supported estimation algorithm.');
    end
    
    if string_match(filter_type, {'kf', 'ekf', 'ukf', 'cdkf', 'srukf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf'})
        process_noise_arg.type           = 'gaussian';
        process_noise_arg.covarianceType = 'full';
        process_noise_arg.tag            = 'GSSM state noise';
        process_noise_arg.dimension      = model.processNoiseDimension;
        process_noise_arg.mean           = zeros(model.processNoiseDimension, 1);
        process_noise_arg.covariance     = [(1e-4*eye(3)).^2 zeros(3); zeros(3) (1e-7*eye(3)).^2]; % [km^2  n/a; n/a (km/sec)^2]
        
        x_noise_infer = generate_noise_model(process_noise_arg);
    end
    
    % Call inference algorithm / estimator
    x_ray_nav_system = X_RayNavigationSystem.create(filter_type, x_nav_inference_model, x_noise_infer, z_noise_infer, time_data_x_ray);
    x_est = x_ray_nav_system.resolve(x_init, z, filter_params);
    
    err_arr(i, 1, :) = (sqrt(sum((x_est(1:3, 2 : end)).^2, 1)) - sqrt(sum((x(1:3, 2 : end)).^2, 1))); % |r| error, ie norm (module) of attitute vector error, [km]
    err_arr(i, 2, :) = 1e3*(sqrt(sum((x_est(4:6, 2 : end)).^2, 1)) - sqrt(sum((x(4:6, 2 : end)).^2, 1))); % |v| error, ie (module) of velocity vector error, [m / sec]
    errors_on_iteration = [squeeze(err_arr(i, 1, :))'; squeeze(err_arr(i, 2, :))';];
    errors_on_iteration_by_coord = 1e3*(x_est(:, 2 : end) - x(:, 2 : end)); % [ [m]; [m]; [m]; [m / sec]; [m / sec]; [m / sec]]
    
    if (draw_iterations)
        figure(1);
        draw_error_on_iteration(time_data_x_ray, errors_on_iteration(1, :), errors_on_iteration(2, :), filter_type, ...
            {'Attitude error, km', 'Velocity error, m/sec'}, {'r'}, {'v'}, 'X-Ray navigation');
        
        figure(2);
        draw_error_on_iteration(time_data_x_ray, errors_on_iteration_by_coord(1:3, :)/1e3, errors_on_iteration_by_coord(4:6, :), ...
            filter_type, {'Attitude error, km', 'Velocity error, m/sec'}, {'r_x', 'r_y', 'r_z'}, {'v_x', 'v_y', 'v_z'}, 'X-Ray navigation');
        
        %{
        z_pred_error = 1e3*(z_est(:, 2:end) - z(:, 2:end));
        figure(3);
        clf;
        plot(time_data_x_ray.TimeInHour(2:end), z_pred_error, 'LineWidth', 1);
        xlabel('time, hours');
        ylabel('z preditcion error, ms');
        hold on;
        grid on;
        title([filter_type ': Observation preditcion error X-Ray navigation']);
        drawnow
            %}
    end
    
    fprintf('iteration %d completed\n', i);
end
toc;

if (number_of_runs > 1)
    calculate_and_drow_estimate_stats(err_arr, time_data_x_ray, {'attitude', 'velocity'}, filter_type, 'X-Ray navigation', {'km', 'm / sec'});
end
disp(' ');
