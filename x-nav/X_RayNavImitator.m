classdef X_RayNavImitator
    properties (SetAccess = private)
        inference_model;
    end
    
    methods
        function obj = X_RayNavImitator(inference_model)
            obj.inference_model = deep_clone(inference_model);
        end
        
        function [state, observations] = resolve(this, initial_state, time_data)
            % random permutation
            rng(sum(100*clock), 'v5normal');
            rng(sum(100*clock), 'v5uniform');
            
            datapoints_count = time_data.SimulationNumber;
            
            state = zeros(this.inference_model.stateDimension, datapoints_count);
            observations = zeros(this.inference_model.observationDimension, datapoints_count);
            
            x_noise = this.inference_model.processNoise.sample(this.inference_model.processNoise, datapoints_count);
            z_noise = this.inference_model.observationNoise.sample(this.inference_model.observationNoise, datapoints_count);
            
            state(:, 1) = initial_state;
            observations(:, 1) = this.inference_model.observation_fun(this.inference_model, state(:, 1), z_noise(1), []);
            
            for j = 2:datapoints_count
                current_time = time_data.Time(j);
                current_epoch = time_data.get_current_epoch(time_data.Time(j));
                
                this.inference_model = this.inference_model.set_params(this.inference_model, current_epoch, current_time);
                state(:, j) = this.inference_model.transition_fun(this.inference_model, state(:, j-1), x_noise(j), []);
                observations(:, j) = this.inference_model.observation_fun(this.inference_model, state(:, j), z_noise(j), []);
            end
        end
    end
    
end
