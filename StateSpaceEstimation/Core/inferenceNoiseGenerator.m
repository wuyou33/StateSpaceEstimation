function [ processNoise, observationNoise, outputInferenceDataStructure ] = inferenceNoiseGenerator( inferenceDataStructure, ...
    estimatorType, ...
    processNoiseAdaptMethod, ...
    processNoiseAdaptParams, ...
    observationNoiseAdaptMethod, ...
    observationNoiseAdaptParams )
%% Generate process and observation noise data structures for a given InferenceDS data structure
% and algorithm type. All ReBEL estimation algorithms take an inference data structure (InferenceDS),
% as well as two system noise data structures (process noise and observation noise) as arguments.
% INPUT
%          InferenceDS         (InferenceDS) Inference data structure generated from a GSSM file by 'geninfds'
%          estimatorType       (string) type of estimator to be used (i.e. 'kf', 'ukf', 'ekf', 'pf', etc.)
%          pNoiseAdaptMethod  <<optional>> (string) Process noise covariance adaptation method :
%                                      'anneal'        : annealing
%                                      'lambda-decay'  : RLS like lambda decay
%                                      'robbins-monro' : Robbins-Monro stochastic approximation
%                               If this field is set, then pNoiseAdaptParams must also be set.
%          pNoiseAdaptParams  <<optional>> (vector) noise adaptation parameters. Depend on pNoiseAdaptMethod
%                                 if 'anneal'        : [annealing_factor minimum_allowed_variance]
%                                 if 'lambda-decay'  : [lambda_factor minimum_allowed_variance]
%                                 if 'robbins-monro' : [1/nu_initial 1/nu_final]
%          oNoiseAdaptMethod  <<optional>> Observation noise covariance adaptation method : same as above
%                                          except the only allowed method is 'robbins-monro'
%          oNoiseAdaptParams  <<optional>> Same as above for process noise
%
%   OUTPUT
%          pNoise              (NoiseDS) process noise data structure
%          oNoise              (NoiseDS) observation noise data structure
%          InferenceDS         (InferenceDS) updated inference data structure

%%
    if ((nargin < 2) || rem(nargin, 2))
        error(' [ inferenceNoiseGenerator ] Not enough input parameters.');
    end

    if (nargout ~= 3)
        error(' [ inferenceNoiseGenerator ] Not enough output arguments.');
    end

    %% 
    outputInferenceDataStructure = deepClone(inferenceDataStructure);

    outputInferenceDataStructure.estimatorType = estimatorType;

    switch (inferenceDataStructure.inferenceType)
        case 'state'               
            processNoise      = deepClone(inferenceDataStructure.model.processNoise);
            observationNoise  = deepClone(inferenceDataStructure.model.observationNoise);

            % Sigma point filters family
            if stringmatch(estimatorType, {'kf','ekf','ukf','cdkf','srukf','srcdkf', 'ckf', 'sckf'})
                % If default noise source is not Guassian, define a Gaussian noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, {'gaussian'})
                    processesNoiseArg.type       = 'gaussian';
                    processesNoiseArg.dimension  = processNoise.dimension;
                    processesNoiseArg.mean       = processNoise.mean;
                    processesNoiseArg.covariance = processNoise.covariance;

                    processNoise = generateNoiseDataSet(processesNoiseArg);
                end

                if ~stringmatch(observationNoise.noiseSourceType, {'gaussian'})
                    observationNoiseArg.type       = 'gaussian';
                    observationNoiseArg.dimension  = observationNoise.dimension;
                    observationNoiseArg.mean       = observationNoise.mean;
                    observationNoiseArg.covariance = observationNoise.covariance;

                    observationNoise = generateNoiseDataSet(observationNoiseArg);
                end

                if stringmatch(estimatorType, {'srukf', 'srcdkf'})
                    processNoise = convertGassianNoise(processNoise, 'sqrt');                
                    observationNoise = convertGassianNoise(observationNoise, 'sqrt');
                end
                if stringmatch(estimatorType, {'sckf'})
                    processNoise = convertGassianNoise(processNoise, 'svd');
                    observationNoise = convertGassianNoise(observationNoise, 'svd');
                end                
            end

            % Particle filter family
            if stringmatch(estimatorType, 'gspf')            
                % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, 'gmm')                
                    processesNoiseArg.type           = 'gmm';
                    processesNoiseArg.covarianceType = 'sqrt'; % GSPF use square-root covariance matrices
                    processesNoiseArg.dimension      = processNoise.dimension;          % process noise dimension                
                    processesNoiseArg.weights        = 1;
                    processesNoiseArg.mean = processNoise.mean;                

                    if isfield(processNoise,'covarianceType')
                        covarianceType = processNoise.covarianceType;
                    else
                        covarianceType='full';
                    end

                    switch covarianceType
                        case {'sqrt','sqrt-diag'}
                            processesNoiseArg.covariance(:,:,1) = processNoise.covariance;
                        case {'full','diag'}
                            processesNoiseArg.covariance(:,:,1) = chol(processNoise.covariance)';
                        otherwise
                            error('[ inferenceNoiseGenerator::gspf ] Unknown process noise covariance type.');
                    end

                    processNoise = generateNoiseDataSet(processesNoiseArg);                
                else                
                    % Make sure the GMM component densities is of cov_type 'sqrt'
                    switch (processNoise.covarianceType)
                        case 'diag'
                            processNoise = convertGassianNoise(processNoise,'sqrt-diag');                        
                        case 'full'
                            processNoise = convertGassianNoise(processNoise,'sqrt');
                    end                
                end
            end

            % 'Gaussian Mixture Sigma-Point Particle Filter'        
            if stringmatch(estimatorType, 'gmsppf')            
                % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, 'gmm')
                    processesNoiseArg.type           = 'gmm';
                    processesNoiseArg.covarianceType = 'sqrt';
                    processesNoiseArg.dimension      = processNoise.dimension;                
                    processesNoiseArg.weights        = 1;
                    processesNoiseArg.mean           = processNoise.mean;

                    if isfield(processNoise,'covarianceType')
                        covarianceType = processNoise.covarianceType;
                    else
                        covarianceType='full';
                    end

                    switch covarianceType
                        case {'sqrt','sqrt-diag'}
                            processesNoiseArg.cov(:,:,1) = processNoise.cov;
                        case {'full','diag'}
                            processesNoiseArg.cov(:,:,1) = chol(processNoise.cov)';
                        otherwise
                            error(' [ inferenceNoiseGenerator::gmsppf ] Unknown process noise covariance type.');
                    end

                    processNoise = generateNoiseDataSet(processesNoiseArg);
                else
                    % Make sure the GMM component densities is of cov_type 'sqrt'
                    switch (processNoise.covarianceType)
                        case 'diag'
                            processNoise = convertGassianNoise(processNoise,'sqrt-diag');
                        case 'full'
                            processNoise = convertGassianNoise(processNoise,'sqrt');
                    end
                end

                % If observation noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance

                if ~stringmatch(observationNoise.noiseSourceType, 'gmm')
                    processesNoiseArg.type           = 'gmm';
                    processesNoiseArg.covarianceType = 'sqrt';
                    processesNoiseArg.dimension      = observationNoise.dimension;
                    processesNoiseArg.weights        = 1;
                    processesNoiseArg.mean           = observationNoise.mean;

                    if isfield(observationNoise,'covarianceType')
                        covarianceType = observationNoise.covarianceType;
                    else
                        covarianceType='full';
                    end

                    switch covarianceType
                        case {'sqrt','sqrt-diag'}
                            processesNoiseArg.cov(:,:,1) = observationNoise.cov;
                        case {'full','diag'}
                            processesNoiseArg.cov(:,:,1) = chol(observationNoise.cov)';
                        otherwise
                            error(' [ inferenceNoiseGenerator::gmsppf ] Unknown observation noise covariance type.');
                    end

                    observationNoise = generateNoiseDataSet(processesNoiseArg);      % generate process noise data structure
                else
                    % Make sure the GMM component densities is of cov_type 'sqrt'
                    %-- process noise
                    switch (observationNoise.generateNoiseDataSet)                        % Determine cov_type of Gaussian noise source
                        case 'diag'
                            observationNoise = convertGassianNoise(observationNoise,'sqrt-diag');
                        case 'full'
                            observationNoise = convertGassianNoise(observationNoise,'sqrt');
                    end
                end

            end        
            processNoise.tag     = 'state';
            observationNoise.tag = 'observation';
        case 'parameter'
            error(' [ inferenceNoiseGenerator] parameter estimation not implelemented');

        case 'joint'
            error(' [ inferenceNoiseGenerator] joint estimation not implelemented');

        otherwise
            error([' [ inferenceDataGenerator ] Inference type ''' args.type ''' not supported.']);
    end

    outputInferenceDataStructure.processNoiseAdaptMethod     = [];
    outputInferenceDataStructure.observationNoiseAdaptMethod = [];
    processNoise.adaptMethod                = [];
    observationNoise.adaptMethod            = [];

    if (nargin >= 4)
        processNoise.adaptMethod = processNoiseAdaptMethod;
        processNoise.adaptParams = processNoiseAdaptParams;
        outputInferenceDataStructure.processNoiseAdaptMethod = processNoiseAdaptMethod;
    end

    if (nargin == 6)
        observationNoise.adaptMethod = observationNoiseAdaptMethod;
        observationNoise.adaptParams = observationNoiseAdaptParams;
        outputInferenceDataStructure.observationNoiseAdaptMethod = observationNoiseAdaptMethod;
    end

end