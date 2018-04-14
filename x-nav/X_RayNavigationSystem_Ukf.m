classdef X_RayNavigationSystem_Ukf < X_RayNavigationSystem
    methods
        function obj = X_RayNavigationSystem_Ukf(inference_model, x_noise_infer, z_noise_infer, time_data)
            obj@X_RayNavigationSystem(inference_model, x_noise_infer, z_noise_infer, time_data);
        end
        
        function state_est = resolve(this, state_initial, z, filter_params)
            alpha = filter_params.alpha; % scale factor (UKF parameter) 1e-3
            beta  = filter_params.beta; % 2 is a optimal setting for Gaussian priors (UKF parameter)
            kappa = filter_params.kappa; % 0 is optimal for state dimension = 2 (UKF parameter)
            
            datapoints_count = this.time_data.SimulationNumber;
            
            this.inference_model.spkfParams = [alpha beta kappa];
            state_est = zeros(this.inference_model.stateDimension, datapoints_count);
            state_est(:, 1) = state_initial;
            x_cov_est = filter_params.x_cov_initial;
            
            for j = 2 : datapoints_count
                current_time = this.time_data.Time(j);
                current_epoch = this.time_data.get_current_epoch(this.time_data.Time(j));
                
                this.inference_model.model = this.inference_model.model.set_params(this.inference_model.model, current_epoch, current_time);
                
                [state_est(:, j), x_cov_est] = ukf(state_est(:, j-1), x_cov_est, this.x_noise_infer, this.z_noise_infer, z(:, j), this.inference_model, [], []);
            end
        end
    end
end
