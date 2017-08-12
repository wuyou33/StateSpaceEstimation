function [ prior, likelihood, evidence, posterior ] = gmmProbability( gmmSet, dataSet, evidenceWeights )
    % gmmProbability. Calculates any of the related (through Bayes rule) probabilities of a Gaussian Mixture Model (gmmSet) and a given dataset (dataSet).
    %       Probabilities are:
    %                   P(X|C) . P(C)                   likelihood . prior
    %       P(C|X) = -----------------   posterior =  --------------------
    %                      P(X)                             evidence
    %
    %   where C is the component classes (Gaussians) of the GMM and X is the data.
    %
    %   [ prior, likelihood, evidence, posterior ] = gmmProbability( gmmSet, dataSet, evidenceWeights )
    %
    %   INPUT:
    %       gmmSet                 Gaussian mixture model data structure with the following fields;
    %           .covarianceType    covariance matrix type ('full' , 'diag' , 'sqrt' , 'sqrt-diag', 'svd');
    %           .dimension         data dimension;
    %           .mixtureCount      number of Gaussian component densities;
    %           .weights           mixing priors (component weights);
    %           .mean              Gaussian component means;
    %           .covariance        covariance matrices of Gaussian components;
    %       dataSet                buffer of N dim-by-1 data set vectors to be evaluated;
    %       evidenceWeights        <optional> vector of sample weights. If specified, the sampleset will be weighted according to these weights.
    %
    %   OUTPUT:
    %       prior          the prior (without seeing data) probability of a componentdensity generating any given data vector, i.e. P(C(i)).
    %                        This is simply the same as the prior mixing weights, 'gmmSet.weights';
    %       likelihood     martrix where the j, i-th entry is the likelihoodof input column vector i (of dataSet) conditioned on component
    %                        density j, i.e. P(X(i)|C(j));
    %       evidence       matrix where the i-th entry is the total data probabilityfor a given data vector X(i), i.e. P(X(i)) = sum_over_all_j [P(X(i)|C(j)) ];
    %       posterior      matrix where the j,i-th entry is the posteriorprobability (after seeing the data) that a component
    %                        density j has generated a specific data vector X(i), i.e. P(C(j)|X(i)) (class posterior probabilities).
    %
    narginchk(2, 3);
    
    [dim, bucketSize] = size(dataSet);
    
    if dim ~= gmmSet.dimension
        error('[ gmmProbability ] Data dimension and GMM model dimension is not the same.');
    end
    
    normFact = (2*pi)^(gmmSet.dimension / 2);
    
    prior = gmmSet.weights(:);
    
    likelihood = zeros(gmmSet.mixtureCount, bucketSize);
    
    for i = 1 : gmmSet.mixtureCount
        centered = dataSet - cvecrep(gmmSet.mean(:, i), bucketSize);
        
        switch gmmSet.covarianceType
            case {'full', 'diag'}
                sqrtCov = chol(gmmSet.covariance(:, :, i), 'lower');
            case {'sqrt', 'sqrt-diag', 'svd'}
                sqrtCov = gmmSet.covariance(:, :, i);
            otherwise
                error(['[ gmmProbability ] Unknown covariance type ', gmmSet.covarianceType]);
        end
        
        x = sqrtCov \ centered;
        likelihood(i, :) = exp(-0.5*sum(x.*x, 1)) / abs( (normFact*prod(diag(sqrtCov))) );
    end
    
    likelihood = likelihood + 1e-99;
    
    if (nargin == 3)
        evidence = prior' * (likelihood ./ evidenceWeights(ones(gmmSet.mixtureCount, 1),:)) + 1e-99;
    else
        evidence = prior' * likelihood + 1e-99;
    end
    
    posterior = likelihood ./ ( (1 ./ prior) * evidence ) + 1e-99;
    posterior = posterior ./ rvecrep(sum(posterior, 1), gmmSet.mixtureCount);
    
%     lll = likelihoodGaussian(gmmSet, dataSet, 0);
end
