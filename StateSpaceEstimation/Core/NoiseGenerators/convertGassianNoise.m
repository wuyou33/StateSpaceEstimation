function noiseDataStructure = convertGassianNoise(noiseDataStructure, targetCovarianceType)
% convertGassianNoise: Convert Gaussian noise source from one covariance type to another

%%

if (nargin < 2)
    error(' [ convertGassianNoise ] Incorrect number of input arguments.');
end

if ~stringmatch(noiseDataStructure.noiseSourceType, {'gaussian','gmm'})
    error(' [ convertGassianNoise ] This function can only operate on Gaussian or Combination-Gaussian noise sources.');
end

if ~ischar(targetCovarianceType)
    error(' [ convertGassianNoise ] Second input argument must be a string.');
end

if ~stringmatch(targetCovarianceType,{'full','diag','sqrt','sqrt-diag'})
    error([' [ convertGassianNoise ] Unknown target cov_type ''' targetCovarianceType '''']);
end

switch targetCovarianceType
    
    case 'full'        
        switch noiseDataStructure.covarianceType
            case {'sqrt','sqrt-diag'}                
                if stringmatch(noiseDataStructure.noiseSourceType, 'gmm')
                    for k=1:noiseDataStructure.mixtureNumber, % M
                        noiseDataStructure.covariance(:,:,k) = noiseDataStructure.covariance(:,:,k) * noiseDataStructure.covariance(:,:,k)';
                    end
                else
                    noiseDataStructure.covariance = noiseDataStructure.covariance * noiseDataStructure.covariance';
                end
        end
        noiseDataStructure.covariance = targetCovarianceType;
        
    case 'diag'        
        % todo: continue here
        switch noiseDataStructure.cov_type
            case {'sqrt','sqrt-diag'}
                if stringmatch(noiseDataStructure.ns_type,'gmm')
                    for k=1:noiseDataStructure.M,
                        noiseDataStructure.cov(:,:,k) = diag(diag(noiseDataStructure.cov(:,:,k)*noiseDataStructure.cov(:,:,k)'));
                    end
                else
                    noiseDataStructure.cov = diag(diag(noiseDataStructure.cov*noiseDataStructure.cov'));
                end
            otherwise
                if stringmatch(noiseDataStructure.ns_type,'gmm')
                    for k=1:noiseDataStructure.M,
                        noiseDataStructure.cov(:,:,k) = diag(diag(noiseDataStructure.cov(:,:,k)));
                    end
                else
                    noiseDataStructure.cov = diag(diag(noiseDataStructure.cov));
                end
        end
        noiseDataStructure.cov_type = targetCovarianceType;
        
        %........................................................................................
    case 'sqrt'
        
        switch noiseDataStructure.cov_type
            case {'full','diag'}
                if stringmatch(noiseDataStructure.ns_type,'gmm')
                    for k=1:noiseDataStructure.M,
                        noiseDataStructure.cov(:,:,k) = chol(noiseDataStructure.cov(:,:,k))';
                    end
                else
                    noiseDataStructure.cov = chol(noiseDataStructure.cov)';
                end
        end
        noiseDataStructure.cov_type = targetCovarianceType;
        
        %........................................................................................
    case 'sqrt-diag'
        
        switch noiseDataStructure.cov_type
            case {'full','diag'}
                if stringmatch(noiseDataStructure.ns_type,'gmm')
                    for k=1:noiseDataStructure.M,
                        noiseDataStructure.cov(:,:,k) = diag(diag(chol(noiseDataStructure.cov(:,:,k))'));
                    end
                else
                    noiseDataStructure.cov = diag(diag(chol(noiseDataStructure.cov)'));
                end
            otherwise
                if stringmatch(noiseDataStructure.ns_type,'gmm')
                    for k=1:noiseDataStructure.M,
                        noiseDataStructure.cov(:,:,k) = diag(diag(noiseDataStructure.cov(:,:,k)));
                    end
                else
                    noiseDataStructure.cov = diag(diag(noiseDataStructure.cov));
                end
        end
        noiseDataStructure.cov_type = targetCovarianceType;
        
        %........................................................................................
    otherwise
        error([' [ convgausns ] Unknown target cov_type ''' targetCovarianceType '''']);
end