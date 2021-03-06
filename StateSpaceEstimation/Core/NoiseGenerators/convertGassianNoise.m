function noiseDataStructure = convertGassianNoise(noiseDataStructure, targetCovarianceType)
    % convertGassianNoise. Convert Gaussian noise source from one covariance type to another (targetCovarianceType).
    %
    % noiseDataStructure = convertGassianNoise(noiseDataStructure, targetCovarianceType)
    %
    %   INPUT
    %       noiseDataStructure      structure which fully describe input noise;
    %       targetCovarianceType    target type (type of covariance matrix) to convert source covariance.
    %
    %   OUTPUT
    %       noiseDataStructure    structure which fully describe noise with converted covariance.
    %
    %% error checking
    if nargin < 2
        error('[ convertGassianNoise ] Incorrect number of input arguments.');
    end
    
    if ~stringmatch(noiseDataStructure.noiseSourceType, {'gaussian', 'gmm'})
        error('[ convertGassianNoise ] This function can only operate on Gaussian or Combination-Gaussian noise sources.');
    end
    
    if ~ischar(targetCovarianceType)
        error('[ convertGassianNoise ] Second input argument must be a string.');
    end
    
    if ~stringmatch(targetCovarianceType, {'full', 'diag', 'sqrt', 'sqrt-diag'})
        error(['[ convertGassianNoise ] Unknown target covarianceType ''' targetCovarianceType '''']);
    end
    %%
    switch targetCovarianceType
        case 'full'
            switch noiseDataStructure.covarianceType
                case {'sqrt', 'sqrt-diag'}
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = noiseDataStructure.covariance(:, :, k) * noiseDataStructure.covariance(:, :, k)';
                        end
                    else
                        noiseDataStructure.covariance = noiseDataStructure.covariance * noiseDataStructure.covariance';
                    end
            end
            noiseDataStructure.covarianceType = targetCovarianceType;
            
        case 'diag'
            switch noiseDataStructure.covarianceType
                case {'sqrt', 'sqrt-diag'}
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = diag(diag(noiseDataStructure.covariance(:, :, k)*noiseDataStructure.covariance(:, :, k)'));
                        end
                    else
                        noiseDataStructure.covariance = diag(diag(noiseDataStructure.covariance*noiseDataStructure.covariance'));
                    end
                otherwise
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = diag(diag(noiseDataStructure.covariance(:, :, k)));
                        end
                    else
                        noiseDataStructure.covariance = diag(diag(noiseDataStructure.covariance));
                    end
            end
            noiseDataStructure.covarianceType = targetCovarianceType;
            
        case 'sqrt'
            switch noiseDataStructure.covarianceType
                case {'full', 'diag'}
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = chol(noiseDataStructure.covariance(:, :, k))';
                        end
                    else
                        if ~isscalar(noiseDataStructure.covariance)
                            noiseDataStructure.covariance = chol(noiseDataStructure.covariance)';
                        end
                    end
            end
            noiseDataStructure.covarianceType = targetCovarianceType;
            
        case 'svd'
            error('svd covariance type does not supported at initialization');
            
        case 'sqrt-diag'
            switch noiseDataStructure.covarianceType
                case {'full', 'diag'}
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = diag(diag(chol(noiseDataStructure.covariance(:, :, k))'));
                        end
                    else
                        noiseDataStructure.covariance = diag(diag(chol(noiseDataStructure.covariance, 'lower')));
                    end
                otherwise
                    if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                        for k = 1:noiseDataStructure.mixtureNumber
                            noiseDataStructure.covariance(:, :, k) = diag(diag(noiseDataStructure.covariance(:, :, k)));
                        end
                    else
                        noiseDataStructure.covariance = diag(diag(noiseDataStructure.covariance));
                    end
            end
            noiseDataStructure.covarianceType = targetCovarianceType;
            
        otherwise
            error(['[ convertGassianNoise ] Unknown target covarianceType ''' targetCovarianceType '''']);
    end
end
