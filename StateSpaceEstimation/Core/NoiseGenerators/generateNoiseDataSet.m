function noiseDataSet = generateNoiseDataSet( args )
    % generateNoiseDataSet. Generates a noiseDataSet data structure describing a noise source (args).
    %
    %   noiseDataSet = generateNoiseDataSet(args)
    %
    %   This function generates noise source which is encapsulated in the noiseDataSet data structure.
    %   All arguments to the function are passed via the args argument data structure which has the
    %   following fields, depending on the required noise source type:
    %
    %   1) Gaussian Noise Source (single Gaussian noise source)
    %       args                   (structure)
    %           .type              (string)   'gaussian';
    %           .covarianceType    (string)   type of covariance;
    %           .tag               (string)   identifier tag ( by default is empty);
    %           .dimension         (scalar)   noise vector length;
    %           .mean              (c-vector) mean vector;
    %           .covariance        (matrix)   covariance matrix (should comply with covarianceType).
    %
    %   2) Scalar Gamma noise source
    %       args                   (structure)
    %           .type              (string)  'gamma';
    %           .tag               (string)  identifier tag;
    %           .covarianceType    (string)  type of covariance;
    %           .dimension         (scalar)  1 (multivariate Gamma noise not yet supported);
    %           .alpha             (scalar)  alpha parameter;
    %           .beta              (scalar)  beta parameter.
    %
    %   The generated noise source data structure, noiseDataSet, has the following required fields.
    %   Depending on the noise source type, the data structure may also have other type dependent fields.
    %
    %     noiseDataSet             (structure)
    %          .noiseSourceType    (string) noise source type;
    %          .dimension          (scalar) noise source dimension;
    %          .sample             (function handle) function to generate N noise source samples;
    %          .likelihood         (function handle) function to evaluate the likelihood of a given noise source sample;
    %          .update             (function handle) <<optional>> function to update the internal structure of noise source;
    %          .covarianceType     (string) type of covariance.
    
    %% error checking
    if nargin ~= 1
        error('[ generateNoiseDataSet ] Incorrect number of inputs');
    end
    
    if ~isstruct(args)
        error('[ generateNoiseDataSet ] The input argument to this function must be an argument data structure.');
    end
    
    if ~(isfield(args, 'type') && ischar(args.type))
        error('[ generateNoiseDataSet ] The argument data structure must have a ''type'' field (string), specifying the desired noise source type.');
    else
        noise.type = args.type;
    end
    
    if isfield(args, 'tag')
        if ischar(args.tag)
            noise.tag = args.tag;
        else
            error('[ generateNoiseDataSet ] The ''tag'' field (must be a string).');
        end
    else
        noise.tag = '';
    end
    
    if (isfield(args, 'dimension') && isnumeric(args.dimension) && (length(args.dimension(:)) == 1 ))
        noise.dimension  = args.dimension;
    else
        error('[ generateNoiseDataSet ] Noise source dimension not specified or not a scalar.');
    end
    
    %% build noise source
    switch (noise.type)
        case 'gamma'
            if noise.dimension ~= 1
                error('[ generateNoiseDataSet::gamma ] Multivariate Gamma noise not supported yet.');
            end
            
            if isfield(args, 'alpha') && isnumeric(args.alpha)
                noise.alpha = args.alpha;
            else
                error('[ generateNoiseDataSet::gamma ] Alpha parameter must be a scalar');
            end
            
            if isfield(args, 'beta') && isnumeric(args.beta)
                noise.beta = args.beta;
            else
                error('[ generateNoiseDataSet::gamma ] Beta parameter must be a scalar');
            end
            
            noiseDataSet.noiseSourceType = noise.type;
            noiseDataSet.dimension       = noise.dimension;
            noiseDataSet.alpha           = noise.alpha;
            noiseDataSet.beta            = noise.beta;
            noiseDataSet.mean            = noise.alpha * noise.beta;
            noiseDataSet.covariance      = noise.alpha * (noise.beta^2);
            noiseDataSet.sample          = @sampleGamma;
            noiseDataSet.likelihood      = @likelihoodGamma;
            noiseDataSet.covarianceType  = args.covarianceType;
            
        case 'gaussian'
            if isfield(args, 'covarianceType')
                if (ischar(args.covarianceType) && stringmatch(args.covarianceType, {'full', 'diag', 'sqrt', 'sqrt-diag'}))
                    noiseDataSet.covarianceType = args.covarianceType;
                else
                    error('[ generateNoiseDataSet::gaussian ] Noise source covarianceType not recognized or not a string.');
                end
            else
                warning('[ covarianceType::gaussian ] Covariance type struct.covarianceType not assigned!. Assuming default value, ''full''');
                noiseDataSet.covarianceType = 'full';
            end
            
            if ~isfield(args, 'mean')
                noise.mean = zeros(noise.dimension, 1);
            else
                noise.mean = args.mean;
            end
            
            if ~isfield(args, 'covariance')
                error('[ generateNoiseDataSet::gaussian ] Covariance struct .covariance not assigned!.');
            end
            
            switch (noiseDataSet.covarianceType)
                case 'full'
                    noiseDataSet.covariance = args.covariance;
                case 'diag'
                    if (args.covariance == diag(diag(args.covariance)))
                        noiseDataSet.covariance = args.covariance;
                    else
                        error(['[ generateNoiseDataSet::gaussian::diag ] Diagonal Guassian noise source cannot have non-zero off diagonal values',...
                            'in the covariance matrix.']);
                    end
                case 'sqrt'
                    if (isnumeric(args.cov) && isequal(size(args.covariance), [noiseDataSet.dimension noiseDataSet.dimension]))
                        noiseDataSet.covariance = args.covariance;
                    else
                        error('[ generateNoiseDataSet::gaussian::sqrt ] Noise source covariance matrix of incorrect dimensions or type.');
                    end
                    
                case 'sqrt-diag'
                    if (isnumeric(args.cov) && isequal(size(args.covariance), [noiseDataSet.dimension noiseDataSet.dimension]))
                        if (args.covariance == diag(diag(args.covariance)))
                            noiseDataSet.covariance = args.covariance;
                        else
                            error(['[ generateNoiseDataSet::gaussian::sqrt-diag ] Diagonal Guassian noise source cannot have non-zero off diagonal values',...
                                'in the covariance matrix.']);
                        end
                    else
                        error('[ generateNoiseDataSet::gaussian::sqrt-diag ] Noise source covariance matrix of incorrect dimensions or type.');
                    end
                otherwise
                    error('[ generateNoiseDataSet::gaussian ] Unknown noise source covarianceType.');
            end
            
            noiseDataSet.noiseSourceType = noise.type;
            noiseDataSet.dimension       = noise.dimension;
            noiseDataSet.mean            = noise.mean;
            noiseDataSet.sample          = @sampleGaussian;
            noiseDataSet.update          = @updateGaussian;
            noiseDataSet.likelihood      = @likelihoodGaussian;
            
        case 'combo-gaussian'
            if (~isfield(args, 'noiseSources') || ~iscell(args.noiseSources))
                error('[ generateNoiseDataSet::combo-gaussian ] Sub noise source field (args.noiseSources) is missing or is not a cell array.');
            end
            
            noise.N = length(args.noiseSources);
            
            if (noise.N < 2)
                error('[ generateNoiseDataSet::combo-gaussian ] A combo-Gaussian noise sources needs at least 2 sub noise sources.');
            end
            
            noiseType = args.noiseSources{1}.noiseSorceType;
            covarianceType = args.noiseSources{1}.covarianceType;
            
            if ~(stringmatch(noiseType,{'gaussian', 'combo-gaussian'}) && stringmatch(covarianceType,{'full', 'diag', 'sqrt', 'sqrt-diag'}))
                error('[ generateNoiseDataSet::combo-gaussian ] A combination Gaussian noise source can only have Gaussian sub noise sources.');
            end
            
            for k = 1:noise.N
                subNoise = args.noiseSources{k};
                if ~stringmatch(subNoise.cov_type,covarianceType)
                    error('[ generateNoiseDataSet::combo-gaussian ] Sub noise sources does not have consistent covTypes: Previous: %s; Current: %s', ...
                        covarianceType, ...
                        subNoise.covarianceType);
                end
            end
            
            noise.covarianceType = covarianceType;
            noise.mean           = zeros(noise.dimension, 1);
            noise.covariance     = zeros(noise.dimension);
            noise.idxArr         = zeros(noise.N, 2);
            
            % Extract sub noise source detail and build combo noise source
            dim = 0;
            for j = 1:noise.N,
                subNoise = args.noiseSources{j};
                dim = dim + subNoise.dimension;
                ind1 = dim - subNoise.dimension + 1;
                ind2 = dim;
                noise.idxArr(j, :) = [ind1 ind2];
                noise.mean(ind1:ind2, 1) = subNoise.mean;
                noise.covariance(ind1:ind2, ind1:ind2) = subNoise.covariance;
            end
            
            if (noise.dimension ~= dimension)
                error('[ generateNoiseDataSet::combo-gaussian ] Combined noise vector dimension does not agree with aggregate dimension of sub noise sources.');
            end
            
            % Restructure noise data structure
            noiseDataSet.type            = 'Noise inference structure';
            noiseDataSet.noiseSourceType = noise.type;
            noiseDataSet.covarianceType  = noise.covarianceType;
            noiseDataSet.tag             = noise.tag;
            noiseDataSet.N               = noise.N;
            noiseDataSet.idxArr          = noise.idxArr;
            noiseDataSet.dimension       = noise.dimension;
            noiseDataSet.mean            = noise.mean;
            noiseDataSet.covariance      = noise.covariance;
            noiseDataSet.noiseSources    = args.noiseSources;
            noiseDataSet.sample          = @sampleComboGaussian;
            noiseDataSet.update          = @updateComboGaussian;
            noiseDataSet.likelihood      = @likelihoodComboGaussian;
            
        case 'combo'
            if (~isfield(args,'noiseSources') || ~iscell(args.noiseSources))
                error('[ generateNoiseDataSet::combo ] Sub noise source field (args.noiseSources) is missing or is not a cell array.');
            end
            
            noise.N = length(args.noiseSources);
            
            if noise.N < 2
                error('[ generateNoiseDataSet::combo ] A combo noise sources needs at least 2 sub noise sources.');
            end
            
            noise.noiseSources = args.noiseSources;
            isGaussian = 1;
            noise.idxArr = zeros(noise.N, 2);
            
            dim = 0;
            for j = 1 : noise.N,
                subNoise = noise.noiseSources{j};
                dim = dim + subNoise.dimension;
                ind1 = dim - subNoise.dimension + 1;
                ind2 = dim;
                noise.idxArr(j,:) = [ind1 ind2];
                
                isGaussian = isGaussian && stringmatch(subNoise.noiseSourceType, {'gaussian', 'combo-gaussian'});
            end
            
            if (noise.dimension ~= dim)
                error('[ generateNoiseDataSet::combo-gaussian ] Combined noise vector dimension does not agree with aggregate dimension of sub noise sources.');
            end
            
            % A Gaussian combination should always be constructed if  the underlying sub noise sources are all Gaussian
            if isGaussian
                arg.type = 'combo-gaussian';
                arg.tag  = noise.tag;
                arg.dimension  = noise.dimension;
                arg.noiseSources = noise.noiseSources;
                noiseDataSet = generateNoiseDataSet(arg);
            else
                noiseDataSet.type            = 'Noise inference structure';
                noiseDataSet.noiseSourceType = noise.type;
                noiseDataSet.covarianceType  = noise.covarianceType;
                noiseDataSet.tag             = noise.tag;
                noiseDataSet.N               = noise.N;
                noiseDataSet.idxArr          = noise.idxArr;
                noiseDataSet.dimension       = noise.dimension;
                noiseDataSet.mean            = noise.mean;
                noiseDataSet.covariance      = noise.covariance;
                noiseDataSet.noiseSources    = noise.noiseSources;
                noiseDataSet.sample          = @sampleCombo;
                noiseDataSet.update          = @updateComboGaussian;
                noiseDataSet.likelihood      = @likelihoodCombo;
            end
            
        case 'gmm' % Gaussian Mixture Model
            if isfield(args, 'mixtureCount')
                noise.mixtureCount = args.mixtureCount;
            else
                error('[ generateNoiseDataSet::gmm ] Number of mixture components not specified.');
            end
            
            if (isequal(size(args.mean), [noise.dimension noise.mixtureCount]))
                noise.mean = args.mean;
            else
                error('[ generateNoiseDataSet::gmm ] Centroid mean dimension error.');
            end
            
            % check for and assign covarianceType
            if isfield(args, 'covarianceType')
                noise.covarianceType = args.covarianceType;
            else
                warning('[ generateNoiseDataSet::gmm ] Covariance type field .covarianceType not assigned!. Assuming default value, ''full''');
                noise.covarianceType = 'full';
            end
            
            if ~isfield(args, 'weights')
                noise.weights = (1 / noise.mixtureCount) * ones(1, noise.mixtureCount);
            else
                if (length(args.weights) == noise.mixtureCount)
                    noise.weights = args.weights / (sum(args.weights));
                else
                    error('[ generateNoiseDataSet::gmm ] Incorrect number of mixing weights (priors).');
                end
            end
            
            switch (noise.covarianceType)
                case {'full', 'diag', 'sqrt', 'sqrt-diag',}
                    if ((noise.mixtureCount == 1) && isequal(size(args.covariance), [noise.dimension noise.dimension]) ||...
                            ((noise.M  > 1) && isequal(size(args.covariance), [noise.dimension noise.dimension noise.mixtureCount])))
                        noise.covariance = args.covariance;
                    else
                        error('[ generateNoiseDataSet::gmm::full ] Noise source covariance matrix has incorrect dimensions.');
                    end
                otherwise
                    error('[ generateNoiseDataSet::gmm ]unknown noise source covarianceType.');
            end
            
            % Restructure noise data structure
            noiseDataSet.type             = 'Noise inference structure';
            noiseDataSet.noiseSourceType  = noise.type;
            noiseDataSet.covarianceType   = noise.covarianceType;
            noiseDataSet.tag              = noise.tag;
            noiseDataSet.dimension        = noise.dimension;
            noiseDataSet.mixtureCount     = noise.mixtureCount;
            noiseDataSet.weights          = noise.weights;
            noiseDataSet.mean             = noise.mean;
            noiseDataSet.covariance       = noise.covariance;
            noiseDataSet.sample           = @gmmSample;
            noiseDataSet.likelihood       = @likelihoodGmm;
            
        otherwise
            error(['[ generateNoiseDataSet ] Noise type ' type ' not supported.']);
    end
    
    noiseDataSet.adaptMethod = [];
end

function noiseDataSet = updateGaussian(noiseDataSet)
    % Updates a Gaussian noise source
    % This function is only a placeholder here since Gaussian noise sources are completely updated, i.e.
    % mean and covariance are set externally, i.e. there are no INTERNAL structure to update.
end
