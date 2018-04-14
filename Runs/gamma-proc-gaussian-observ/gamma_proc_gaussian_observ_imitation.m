clc; close all; clearvars; clear memoize;
addpath(genpath('./'));

fprintf('\nThe demonstration of state space esitmation problem of scalar nonlinear system with Gamma process and Gaussian observation noise.\n');

% {'ukf', 'pf', 'gspf', 'sppf', 'gmsppf'}
filter_types = {'ukf'};

number_of_runs = 40;
datapoints_count  = 50;
draw_iterations = number_of_runs <= 50;
err_arr = zeros(number_of_runs, datapoints_count - 1);

for filter_index = 1 : length(filter_types)
    filter_type = filter_types{filter_index};
    
    disp(filter_type);
    
    model = gssm_gamma_proc_gauss_observ('init');
    
    arg.type    = 'state';
    arg.tag     = 'State estimation for GSSM_N1 system.';
    arg.model   = model;
    arg.algorithm = filter_type;
    
    inference_model = inference_model_generator(arg);
    [proc_noise_infer, obs_noise_infer, inference_model] = inference_noise_generator(inference_model, filter_type);
    
    tic;
    for i = 1 : number_of_runs
        rng(sum(100*clock), 'v5normal');
        rng(sum(100*clock), 'v5uniform');
        
        x  = zeros(model.stateDimension, datapoints_count); % state data buffer
        z  = zeros(model.observationDimension, datapoints_count);   % observation data buffer
        
        x_noise = model.processNoise.sample(model.processNoise, datapoints_count);
        z_noise = model.observationNoise.sample(model.observationNoise, datapoints_count);
        
        z(1) = model.observation_fun(model, x(1), z_noise(1), 1);
        for j = 2:datapoints_count
            x(j) = model.transition_fun(model, x(:,j-1), x_noise(j-1), j-1);
            z(j) = model.observation_fun(model, x(:,j), z_noise(j), j);
        end
        
        u1 = 0 : datapoints_count-1;
        u2 = 1 : datapoints_count;
        
        x_est = zeros(inference_model.stateDimension, datapoints_count);
        x_est(:, 1) = 1; % 21
        x_cov_est = 3 / 4 * eye(inference_model.stateDimension);
        
        % Call inference algorithm / estimator
        switch filter_type
            case 'ukf'
                alpha = 1; % scale factor (UKF parameter)
                beta  = 2; % optimal setting for Gaussian priors (UKF parameter)
                kappa = 0; % optimal for state dimension=2 (UKF parameter)
                
                inference_model.spkfParams = [alpha beta kappa];
                
                for k = 2 : datapoints_count
                    [x_est(:, k), x_cov_est] = ukf(x_est(:, k-1), x_cov_est, proc_noise_infer, obs_noise_infer, z(k), inference_model, u1(k-1), u2(k));
                end
            case 'pf'
                particlesCount = 5e4; % 200
                particle_set.particlesNum = particlesCount;
                particle_set.particles = randn(inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
                particle_set.weights = column_vector_replicate(1 / particlesCount, particlesCount);
                inference_model.resampleThreshold = 0.1;
                inference_model.estimateType = 'mean';
                inference_model.resampleMethod = 'residual-kitagawa';
                
                for k = 2 : datapoints_count
                    [x_est(:, k), particle_set] = pf(particle_set, proc_noise_infer, obs_noise_infer, z(k), u1(k-1), u2(k), inference_model);
                end
            case 'gspf'
                particlesCount = 10000; % 200
                particle_set.particlesNum = particlesCount;
                initial_particles = randn(inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
                
                % fit a n component GMM to initial state distribution
                particle_set.stateGMM = gmm_fit(initial_particles, 4, [0.001 10], 'sqrt');
                
                inference_model.estimateType = 'mean';
                inference_model.threshold = 0.001;
                inference_model.resampleMethod = 'residual';
                
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
                
                for k = 2 : datapoints_count
                    [x_est(:, k), particle_set] = gspf(particle_set, gspf_process_noise, obs_noise_infer, z(k), u1(k-1), u2(k), inference_model);
                end
            case 'sppf'
                particlesCount = 400; % 200
                particle_set.particlesNum = particlesCount;
                particle_set.particles    = column_vector_replicate(x_est(:, 1), particlesCount);
                particle_set.particlesCov = repmat(eye(inference_model.stateDimension), [1 1 particlesCount]);
                
                proc_noise_gauss.covariance = sqrt(2*3/4);
                observ_noise_gauss.covariance = sqrt(1e-1);
                
                [proc_noise_gauss, observ_noise_gauss, ~] = inference_noise_generator(inference_model, 'srukf');
                
                particle_set.processNoise = proc_noise_gauss;
                particle_set.observationNoise = observ_noise_gauss;
                particle_set.weights = column_vector_replicate(1 / particlesCount, particlesCount);
                
                inference_model.spkfType = 'srukf';
                inference_model.spkfParams = [1 2 0]; % [1 0 2];
                inference_model.resampleThreshold = 1;
                inference_model.estimateType = 'mean';
                inference_model.resampleMethod = 'residual';
                
                for k = 2 : datapoints_count
                    [x_est(:, k), particle_set] = sppf(particle_set, proc_noise_infer, obs_noise_infer, z(k), u1(k-1), u2(k), inference_model);
                end
            case 'gmsppf'
                particlesCount = 10000; % 200
                particle_set.particlesNum = particlesCount;
                
                initial_particles = randn(inference_model.stateDimension, particlesCount) + column_vector_replicate(x_est(:, 1), particlesCount);
                
                tempCov = zeros(1, 1, 2);
                tempCov(:, :, 1) = sqrt(2);
                tempCov(:, :, 2) = 1;
                
                particle_set.stateGMM = gmm_fit(initial_particles, 4, [0.001 10], 'sqrt');
                
                inference_model.estimateType = 'mean';
                inference_model.spkfType = 'srukf';
                inference_model.spkfParams = [1 2 0]; % [1 0 2];
                inference_model.resampleMethod = 'residual';
                
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
                
                for k = 2 : datapoints_count
                    [x_est(:, k), particle_set] = gmsppf(particle_set, gmsppf_process_noise, obs_noise_infer, z(k), u1(k-1), u2(k), inference_model);
                end
            otherwise
                error('[ filter_type ] Unknown / not supported estimation algorithm.');
        end
        
        if (draw_iterations)
            figure(1);
            clf;
            p1 = plot(x(1, :), 'LineWidth', 1);
            hold on;
            grid on;
            p2 = plot(z, 'g+', 'LineWidth', 1);
            p3 = plot(x_est(1,:), 'r', 'LineWidth', 1); hold on;
            legend([p1 p2 p3], 'clean', 'noisy', [filter_type ' estimate']);
            xlabel('time');
            title([filter_type ': Nonlinear Time Variant State Estimation (non Gaussian noise)']);
            drawnow
        end
        
        err_arr(i, :) = x_est(1, 2 : end) - x(1, 2:end);
    end
    toc;
    
    mean_x_err = mean(err_arr);
    std_x_err = std(err_arr);
    % sqrt(mean((Xh(1,2:end)-X(1,2:end)).^2))
    rmse_x_err = sqrt( mean( (err_arr).^2 ) );
    figure();
    clf;
    p1 = plot(mean_x_err, 'b', 'LineWidth', 1); hold on;
    p2 = plot(mean_x_err - std_x_err,'r', 'LineWidth', 1); hold on;
    p3 = plot(mean_x_err + std_x_err,'r', 'LineWidth', 1); hold off;
    legend([p1 p2 p3], 'mean rmse','- 1 sigma', '+ 1 sigma');
    grid on;
    xlabel('time');
    title([filter_type ': Nonlinear Time Variant State Estimation (non Gaussian noise)']);
    
    figure();
    clf;
    plot(rmse_x_err, 'LineWidth', 1);
    grid on;
    xlabel('time');
    title([filter_type ': RMSE (non Gaussian noise)']);
    
    disp(' ');
end

