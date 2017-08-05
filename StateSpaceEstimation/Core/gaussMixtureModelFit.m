function [ gmmSet ] = gaussMixtureModelFit( data, mixture, terminationThreshold, covarianceType, arbitraryWidth, evidenceWeights, visualize )
    % gaussMixtureModelFit. Fit a Gaussian mixture model (GMM) with M components to dataset (gmmSet)
    % using an expectation maximization (EM) algorithm.
    %
    %   For more details see
    %       http://stats.stackexchange.com/questions/72774/numerical-example-to-understand-expectation-maximization
    %       http://stackoverflow.com/questions/11808074/what-is-an-intuitive-explanation-of-expectation-maximization-technique
    %       http://stackoverflow.com/questions/15513792/expectation-maximization-coin-toss-examples
    %
    %   [ gmmSet ] = gaussMixtureModelFit( data, mixtureCount, terminationThreshold, covarianceType, evidenceWeights )
    %
    %   INPUT
    %       data                 	Dataset of N samples (column vectors);
    %       mixture              	Number of Gaussian mixture component densities OR a pre-initialized GMM data structure;
    %       terminationThreshold 	Termination threshold 0 < tt < 1 (if change in log likelihood falls below this value, the EM algorithm terminates.)
    %                                   OR if tt is a [1-by-2 vector] the first component is the termination threshold and the second component is the maximum
    %                                   number of iterations allowed for the EM, i.e. thresholdtt iterations];
    %       covarianceType       	Covariance type 'full', 'diag', 'sqrt', 'sqrt-diag', 'svd';
    %       arbitraryWidth       	<optional>  Arbitrary width used if variance collapses to zero: make it 'large'
    %                                          so that centre is responsible for a reasonable number of points;
    %       evidenceWeights      	<optional> vector (1xN) of sample weights used to do a weighted EM fit to the data;
    %       visualize               <optional> boolean value, if true then the function will log error, otherwise will not.
    %
    %   OUTPUT
    %       gmmSet              - Gaussian mixture model data structure with the following fields;
    %           .covarianceType - Covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag';
    %           .dimension      - Data dimension;
    %           .mixtureCount   - Number of Gaussian component densities;
    %           .weights        - Mixing priors (component weights);
    %           .mean           - Gaussian component means;
    %           .covariance     - Covariance matrices of Gaussian components (must comply with cov_type).
    %
    narginchk(4, 7);
    
    [dim, dataDim] = size(data);
    
    % number of EM iterations
    if length(terminationThreshold) == 2
        iterationCount = terminationThreshold(2);
    else
        iterationCount = 100;
    end
    
    if isnumeric(mixture)
        gmmSet.covarianceType = covarianceType;
        gmmSet.dimension = dim;
        gmmSet.mixtureCount = mixture;
        gmmSet.weights = ones(dim, mixture);
        gmmSet.mean = zeros(dim, mixture);
        gmmSet.covariance = repmat(eye(dim), [1 1 mixture]);
    else
        gmmSet = mixture;
        mixture = gmmSet.mixtureCount;
    end
    
    if nargin >= 5
        gmmSet = gaussMixtureModelInitialize(gmmSet, data, iterationCount, arbitraryWidth);
    else
        gmmSet = gaussMixtureModelInitialize(gmmSet, data, iterationCount);
    end
    
    initCovariance = gmmSet.covariance;
    
    prevErr = -Inf;
    
    for i = 1 : iterationCount
        if nargin >= 6
            [~, ~, evidence, posterior] = gmmProbability(gmmSet, data, evidenceWeights);
        else
            [~, ~, evidence, posterior] = gmmProbability(gmmSet, data);
        end
        
        err = sum(log(evidence));
        
        if (nargin == 7 && visualize)
            % error value is negative log likelihood of data evidence
            fprintf(1, '[gaussMixtureModelFit:: Cycle %4d  Evidence %11.6f\n]', n, err);
        end
        
        if (i > 1 && abs((err - prevErr) / prevErr) < terminationThreshold(1))
            return;
        else
            prevErr = err;
        end
        
        paramEst = sum(posterior, 2)';
        
        % Now move new estimates to old parameter vectors
        gmmSet.weights = paramEst / dataDim;
        mean = (posterior * data')' ./ rvecrep(paramEst, dim);
        gmmSet.mean = mean;
        
        for j = 1 : mixture
            centered = data - cvecrep(mean(:, j), dataDim);
            
            switch covarianceType
                case 'full'
                    centered = centered .* rvecrep(sqrt(posterior(j, :)), dim);
                    cov = (centered*centered') / paramEst(j);
                    
                    if min(svd(cov)) < eps
                        cov = initCovariance(:, :, j);
                    end
                case {'sqrt', 'sqrt-diag', 'svd'}
                    centered = centered .* rvecrep((1/sqrt(paramEst(j))) * sqrt(posterior(j, :)), dim);
                    
                    [~, cov] = qr(centered', 0);
                    cov = cov';
                    
                    if min(abs(diag(cov))) < sqrt(eps)
                        cov = initCovariance(:, :, j);
                    end
                case 'diag'
                    dd = sum((centered.*centered).* (rvecrep(posterior(j, :), dim)) , 2) / paramEst(j);
                    cov = diag(dd);
                    
                    if min(dd) < eps
                        cov = initCovariance(:, :, j);
                    end
                otherwise
                    error(['[ gaussMixtureModelFit ] Unknown covariance type ', covarianceType]);
            end
            
            gmmSet.covariance(:, :, j) = cov;
        end
    end
end