function [ state_noise, observation_noise, out_inference_model ] = inference_noise_generator( inference_model, ...
        filter_type, ...
        state_noise_adapt_method, ...
        state_noise_adapt_params, ...
        observation_noise_adapt_method, ...
        observation_noise_adapt_params )
    % inference_noise_generator. Generate process and observation noise data structures for a given inference_model data structure
    % and algorithm type. All estimation algorithms take an inference data structure (inference_model),
    % as well as two system noise data structures (process noise and observation noise) as arguments.
    %
    % [ state_noise, observation_noise, out_inference_model ]
    %   = inference_noise_generator( inference_model, ...
    %       filter_type, ...
    %       state_noise_adapt_method, ...
    %       state_noise_adapt_params, ...
    %       observation_noise_adapt_method, ...
    %       observation_noise_adapt_params )
    %
    %   INPUT
    %       inference_model      inference data structure generated from;
    %       filter_type               type of estimator to be used (i.e. 'kf', 'ukf', 'ekf', 'pf', etc.);
    %       state_noise_adapt_method     <<optional>> Process noise covariance adaptation method (if this field is set, then pNoiseAdaptParams required):
    %                                     'anneal'        : annealing;
    %                                     'lambda-decay'  : RLS like lambda decay;
    %                                     'robbins-monro' : Robbins-Monro stochastic approximation;
    %                                     'sage-husa'     : Sage-Husa estimation;
    %       state_noise_adapt_params     <<optional>> (vector) noise adaptation parameters. Depend on pNoiseAdaptMethod
    %                                      if 'anneal'        : [annealing_factor minimum_allowed_variance]
    %                                      if 'lambda-decay'  : [lambda_factor minimum_allowed_variance]
    %                                      if 'robbins-monro' : [1/nu_initial 1/nu_final];
    %                                      if 'sage-husa'     : sage-Husa estimation method;
    %       observation_noise_adapt_method  <<optional>> observation noise covariance adaptation method : same as process;
    %                                      except the only allowed method is 'robbins-monro';
    %       observation_noise_adapt_params  <<optional>> same as above for process noise.
    %
    %   OUTPUT
    %          state_noise                  - process noise data structure;
    %          observation_noise              - observation noise data structure;
    %          out_inference_model  - updated inference data structure.
    %
    if ((nargin < 2) || rem(nargin, 2))
        error('[ inference_noise_generator ] Not enough input parameters.');
    end
    
    if nargout ~= 3
        error('[ inference_noise_generator ] Not enough output arguments.');
    end
    
    out_inference_model = deep_clone(inference_model);
    
    out_inference_model.filter_type = filter_type;
    
    switch (inference_model.inferenceType)
        case 'state'
            state_noise      = deep_clone(inference_model.model.state_noise);
            observation_noise  = deep_clone(inference_model.model.observation_noise);
            
            % sigma point filters family
            if string_match(filter_type, {'kf', 'ekf', 'ukf', 'cdkf', 'srukf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf'})
                % If default noise source is not Guassian, define a Gaussian noise source with the same dimension, mean and covariance
                if ~string_match(state_noise.noiseSourceType, {'gaussian'})
                    processNoiseArg.type       = 'gaussian';
                    processNoiseArg.dimension  = state_noise.dimension;
                    processNoiseArg.mean       = state_noise.mean;
                    processNoiseArg.covariance = state_noise.covariance;
                    processNoiseArg.covarianceType = 'full';
                    
                    state_noise = generate_noise_model(processNoiseArg);
                end
                
                if ~string_match(observation_noise.noiseSourceType, {'gaussian'})
                    observNoiseArg.type       = 'gaussian';
                    observNoiseArg.dimension  = observation_noise.dimension;
                    observNoiseArg.mean       = observation_noise.mean;
                    observNoiseArg.covariance = observation_noise.covariance;
                    
                    observation_noise = generate_noise_model(observNoiseArg);
                end
                
                if string_match(filter_type, {'srukf', 'srcdkf', 'sckf', 'sppf'})
                    state_noise = convert_gassian_noise(state_noise, 'sqrt');
                    observation_noise = convert_gassian_noise(observation_noise, 'sqrt');
                end
            end
            
            % particle filter family
            if string_match(filter_type, 'gspf')
                % if process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~string_match(state_noise.noiseSourceType, 'gmm')
                    processNoiseArg.type            = 'gmm';
                    processNoiseArg.covarianceType  = 'sqrt';
                    processNoiseArg.dimension       = state_noise.dimension;
                    processNoiseArg.weights         = 1;
                    processNoiseArg.mean            = state_noise.mean;
                    processNoiseArg.mixtureCount    = 1;
                    
                    if isfield(state_noise, 'covarianceType')
                        covarianceType = state_noise.covarianceType;
                    else
                        covarianceType = 'full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag'}
                            processNoiseArg.covariance(:, :, 1) = state_noise.covariance;
                        case {'full', 'diag'}
                            processNoiseArg.covariance(:, :, 1) = chol(state_noise.covariance, 'lower');
                        otherwise
                            error('[ inference_noise_generator::gspf ] Unknown process noise covariance type.');
                    end
                    
                    state_noise = generate_noise_model(processNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (state_noise.covarianceType)
                        case 'diag'
                            state_noise = convert_gassian_noise(state_noise, 'sqrt-diag');
                        case 'full'
                            state_noise = convert_gassian_noise(state_noise, 'sqrt');
                    end
                end
            end
            
            % 'Gaussian Mixture Sigma-Point Particle Filter'
            if string_match(filter_type, 'gmsppf')
                % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~string_match(state_noise.noiseSourceType, 'gmm')
                    processNoiseArg.type           = 'gmm';
                    processNoiseArg.covarianceType = 'sqrt';
                    processNoiseArg.dimension      = state_noise.dimension;
                    processNoiseArg.weights        = 1;
                    processNoiseArg.mean           = state_noise.mean;
                    processNoiseArg.mixtureCount   = 1;
                    
                    if isfield(state_noise, 'covarianceType')
                        covarianceType = state_noise.covarianceType;
                    else
                        warning('[ inference_noise_generator ] covarianceType not defined. meant full by default')
                        covarianceType = 'full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag'}
                            processNoiseArg.covariance(:, :, 1) = state_noise.covariance;
                        case {'full', 'diag'}
                            processNoiseArg.covariance(:, :, 1) = chol(state_noise.covariance, 'lower');
                        otherwise
                            error('[ inference_noise_generator::gmsppf ] Unknown process noise covariance type.');
                    end
                    
                    state_noise = generate_noise_model(processNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (state_noise.covarianceType)
                        case 'diag'
                            state_noise = convert_gassian_noise(state_noise, 'sqrt-diag');
                        case 'full'
                            state_noise = convert_gassian_noise(state_noise, 'sqrt');
                    end
                end
                
                % If observation noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~string_match(observation_noise.noiseSourceType, 'gmm')
                    observNoiseArg.type           = 'gmm';
                    observNoiseArg.covarianceType = 'sqrt';
                    observNoiseArg.dimension      = observation_noise.dimension;
                    observNoiseArg.weights        = 1;
                    observNoiseArg.mean           = observation_noise.mean;
                    observNoiseArg.mixtureCount   = 1;
                    
                    if isfield(observation_noise,'covarianceType')
                        covarianceType = observation_noise.covarianceType;
                    else
                        covarianceType='full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag'}
                            observNoiseArg.covariance(:, :, 1) = observation_noise.covariance;
                        case {'full', 'diag'}
                            observNoiseArg.covariance(:, :, 1) = chol(observation_noise.covariance, 'lower');
                        otherwise
                            error('[ inference_noise_generator::gmsppf ] Unknown observation noise covariance type.');
                    end
                    
                    observation_noise = generate_noise_model(observNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (observation_noise.generate_noise_model)
                        case 'diag'
                            observation_noise = convert_gassian_noise(observation_noise, 'sqrt-diag');
                        case 'full'
                            observation_noise = convert_gassian_noise(observation_noise, 'sqrt');
                    end
                end
                
            end
            
            % 'Sigma point particle filter'
            if (string_match(filter_type, 'sppf'))
                % nothing to change here
            end
            
            state_noise.tag     = 'state';
            observation_noise.tag = 'observation';
            
        case 'parameter'
            error('[ inference_noise_generator::args::type] parameter estimation not implelemented');
            
        case 'joint'
            error('[ inference_noise_generator::args::type] joint estimation not implelemented');
            
        otherwise
            error(['[ inference_model_generator::args::type ] Estimation type ''' args.type ''' not supported.']);
    end
    
    out_inference_model.state_noise_adapt_method     = [];
    out_inference_model.observation_noise_adapt_method = [];
    state_noise.adaptMethod                = [];
    observation_noise.adaptMethod            = [];
    
    if nargin >= 4
        state_noise.adaptMethod = state_noise_adapt_method;
        state_noise.adaptParams = state_noise_adapt_params;
        out_inference_model.state_noise_adapt_method = state_noise_adapt_method;
    end
    
    if (nargin == 6)
        observation_noise.adaptMethod = observation_noise_adapt_method;
        observation_noise.adaptParams = observation_noise_adapt_params;
        out_inference_model.observation_noise_adapt_method = observation_noise_adapt_method;
    end
    
    out_inference_model.model.state_noise = state_noise;
    inference_model.model.observation_noise   = observation_noise;
end
