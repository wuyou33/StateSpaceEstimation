classdef (Abstract) X_RayNavigationSystem < handle
    properties (SetAccess = protected)
        inference_model;
        x_noise_infer;
        z_noise_infer;
        time_data;
    end
    
    methods
        function obj = X_RayNavigationSystem(inference_model, x_noise_infer, z_noise_infer, time_data)
            obj.inference_model = inference_model;
            obj.x_noise_infer = x_noise_infer;
            obj.z_noise_infer = z_noise_infer;
            obj.time_data = time_data;
        end
    end
    
    methods (Abstract)
        state_est = resolve(this, state_initial, z, filter_params);
    end
    
    methods (Static)
        function obj = create(filter_type, inference_model, x_noise_infer, z_noise_infer, time_data)
            switch filter_type
                case 'ukf'
                    obj = X_RayNavigationSystem_Ukf(inference_model, x_noise_infer, z_noise_infer, time_data);
                otherwise
                    error(['[ X_RayNavigationSystem::create ] Inference type ''' filter_type ''' not supported.']);
            end
        end
    end
end
