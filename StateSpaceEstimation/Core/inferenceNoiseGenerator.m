function [ processNoise, observationNoise, outputInferenceDataStructure ] = inferenceNoiseGenerator( inferenceDataStructure, ...
        estimatorType, ...
        processNoiseAdaptMethod, ...
        processNoiseAdaptParams, ...
        observationNoiseAdaptMethod, ...
        observationNoiseAdaptParams )
    % inferenceNoiseGenerator. Generate process and observation noise data structures for a given inferenceDataStructure data structure
    % and algorithm type. All estimation algorithms take an inference data structure (inferenceDataStructure),
    % as well as two system noise data structures (process noise and observation noise) as arguments.
    %
    % [ processNoise, observationNoise, outputInferenceDataStructure ]
    %   = inferenceNoiseGenerator( inferenceDataStructure, ...
    %       estimatorType, ...
    %       processNoiseAdaptMethod, ...
    %       processNoiseAdaptParams, ...
    %       observationNoiseAdaptMethod, ...
    %       observationNoiseAdaptParams )
    %
    %   INPUT
    %       inferenceDataStructure      inference data structure generated from;
    %       estimatorType               type of estimator to be used (i.e. 'kf', 'ukf', 'ekf', 'pf', etc.);
    %       processNoiseAdaptMethod     <<optional>> Process noise covariance adaptation method (if this field is set, then pNoiseAdaptParams required):
    %                                     'anneal'        : annealing;
    %                                     'lambda-decay'  : RLS like lambda decay;
    %                                     'robbins-monro' : Robbins-Monro stochastic approximation;
    %                                     'sage-husa'     : Sage-Husa estimation;
    %       processNoiseAdaptParams     <<optional>> (vector) noise adaptation parameters. Depend on pNoiseAdaptMethod
    %                                      if 'anneal'        : [annealing_factor minimum_allowed_variance]
    %                                      if 'lambda-decay'  : [lambda_factor minimum_allowed_variance]
    %                                      if 'robbins-monro' : [1/nu_initial 1/nu_final];
    %                                      if 'sage-husa'     : sage-Husa estimation method;
    %       observationNoiseAdaptMethod  <<optional>> observation noise covariance adaptation method : same as process;
    %                                      except the only allowed method is 'robbins-monro';
    %       observationNoiseAdaptParams  <<optional>> same as above for process noise.
    %
    %   OUTPUT
    %          processNoise                  - process noise data structure;
    %          observationNoise              - observation noise data structure;
    %          outputInferenceDataStructure  - updated inference data structure.
    %
    if ((nargin < 2) || rem(nargin, 2))
        error('[ inferenceNoiseGenerator ] Not enough input parameters.');
    end
    
    if nargout ~= 3
        error('[ inferenceNoiseGenerator ] Not enough output arguments.');
    end
    
    outputInferenceDataStructure = deepClone(inferenceDataStructure);
    
    outputInferenceDataStructure.estimatorType = estimatorType;
    
    switch (inferenceDataStructure.inferenceType)
        case 'state'
            processNoise      = deepClone(inferenceDataStructure.model.processNoise);
            observationNoise  = deepClone(inferenceDataStructure.model.observationNoise);
            
            % sigma point filters family
            if stringmatch(estimatorType, {'kf', 'ekf', 'ukf', 'cdkf', 'srukf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf'})
                % If default noise source is not Guassian, define a Gaussian noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, {'gaussian'})
                    processNoiseArg.type       = 'gaussian';
                    processNoiseArg.dimension  = processNoise.dimension;
                    processNoiseArg.mean       = processNoise.mean;
                    processNoiseArg.covariance = processNoise.covariance;
                    
                    processNoise = generateNoiseDataSet(processNoiseArg);
                end
                
                if ~stringmatch(observationNoise.noiseSourceType, {'gaussian'})
                    observNoiseArg.type       = 'gaussian';
                    observNoiseArg.dimension  = observationNoise.dimension;
                    observNoiseArg.mean       = observationNoise.mean;
                    observNoiseArg.covariance = observationNoise.covariance;
                    
                    observationNoise = generateNoiseDataSet(observNoiseArg);
                end
                
                if stringmatch(estimatorType, {'srukf', 'srcdkf'})
                    processNoise = convertGassianNoise(processNoise, 'sqrt');
                    observationNoise = convertGassianNoise(observationNoise, 'sqrt');
                end
                if stringmatch(estimatorType, {'sckf', 'fdckf'})
                    processNoise = convertGassianNoise(processNoise, 'svd');
                    observationNoise = convertGassianNoise(observationNoise, 'svd');
                end
            end
            
            % particle filter family
            if stringmatch(estimatorType, 'gspf')
                % if process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, 'gmm')
                    processNoiseArg.type            = 'gmm';
                    processNoiseArg.covarianceType  = 'sqrt';
                    processNoiseArg.dimension       = processNoise.dimension;
                    processNoiseArg.weights         = 1;
                    processNoiseArg.mean            = processNoise.mean;
                    processNoiseArg.mixtureCount    = 1;
                    
                    if isfield(processNoise, 'covarianceType')
                        covarianceType = processNoise.covarianceType;
                    else
                        covarianceType = 'full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag', 'svd'}
                            processNoiseArg.covariance(:, :, 1) = processNoise.covariance;
                        case {'full', 'diag'}
                            processNoiseArg.covariance(:, :, 1) = chol(processNoise.covariance, 'lower');
                        otherwise
                            error('[ inferenceNoiseGenerator::gspf ] Unknown process noise covariance type.');
                    end
                    
                    processNoise = generateNoiseDataSet(processNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (processNoise.covarianceType)
                        case 'diag'
                            processNoise = convertGassianNoise(processNoise, 'sqrt-diag');
                        case 'full'
                            processNoise = convertGassianNoise(processNoise, 'sqrt');
                    end
                end
            end
            
            % 'Gaussian Mixture Sigma-Point Particle Filter'
            if stringmatch(estimatorType, 'gmsppf')
                % If process noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~stringmatch(processNoise.noiseSourceType, 'gmm')
                    processNoiseArg.type           = 'gmm';
                    processNoiseArg.covarianceType = 'sqrt';
                    processNoiseArg.dimension      = processNoise.dimension;
                    processNoiseArg.weights        = 1;
                    processNoiseArg.mean           = processNoise.mean;
                    processNoiseArg.mixtureCount   = 1;
                    
                    if isfield(processNoise, 'covarianceType')
                        covarianceType = processNoise.covarianceType;
                    else
                        warning('[ inferenceNoiseGenerator ] covarianceType not defined. meant full by default')
                        covarianceType = 'full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag', 'svd'}
                            processNoiseArg.covariance(:, :, 1) = processNoise.covariance;
                        case {'full', 'diag'}
                            processNoiseArg.covariance(:, :, 1) = chol(processNoise.covariance, 'lower');
                        otherwise
                            error('[ inferenceNoiseGenerator::gmsppf ] Unknown process noise covariance type.');
                    end
                    
                    processNoise = generateNoiseDataSet(processNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (processNoise.covarianceType)
                        case 'diag'
                            processNoise = convertGassianNoise(processNoise, 'sqrt-diag');
                        case 'full'
                            processNoise = convertGassianNoise(processNoise, 'sqrt');
                    end
                end
                
                % If observation noise source is not a GMM, define a GMM noise source with the same dimension, mean and covariance
                if ~stringmatch(observationNoise.noiseSourceType, 'gmm')
                    observNoiseArg.type           = 'gmm';
                    observNoiseArg.covarianceType = 'sqrt';
                    observNoiseArg.dimension      = observationNoise.dimension;
                    observNoiseArg.weights        = 1;
                    observNoiseArg.mean           = observationNoise.mean;
                    observNoiseArg.mixtureCount   = 1;
                    
                    if isfield(observationNoise,'covarianceType')
                        covarianceType = observationNoise.covarianceType;
                    else
                        covarianceType='full';
                    end
                    
                    switch covarianceType
                        case {'sqrt', 'sqrt-diag', 'svd'}
                            observNoiseArg.covariance(:, :, 1) = observationNoise.covariance;
                        case {'full', 'diag'}
                            observNoiseArg.covariance(:, :, 1) = chol(observationNoise.covariance, 'lower');
                        otherwise
                            error('[ inferenceNoiseGenerator::gmsppf ] Unknown observation noise covariance type.');
                    end
                    
                    observationNoise = generateNoiseDataSet(observNoiseArg);
                else
                    % Make sure the GMM component densities is of covarianceType 'sqrt'
                    switch (observationNoise.generateNoiseDataSet)
                        case 'diag'
                            observationNoise = convertGassianNoise(observationNoise, 'sqrt-diag');
                        case 'full'
                            observationNoise = convertGassianNoise(observationNoise, 'sqrt');
                    end
                end
                
            end
            processNoise.tag     = 'state';
            observationNoise.tag = 'observation';
        case 'parameter'
            error('[ inferenceNoiseGenerator::args::type] parameter estimation not implelemented');
            
        case 'joint'
            error('[ inferenceNoiseGenerator::args::type] joint estimation not implelemented');
            
        otherwise
            error(['[ inferenceDataGenerator::args::type ] Estimation type ''' args.type ''' not supported.']);
    end
    
    outputInferenceDataStructure.processNoiseAdaptMethod     = [];
    outputInferenceDataStructure.observationNoiseAdaptMethod = [];
    processNoise.adaptMethod                = [];
    observationNoise.adaptMethod            = [];
    
    if nargin >= 4
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
