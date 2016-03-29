function [prior, likelihood, evidence, posterior] = gmmProbability( gmmDataSet, x, w )
% GMMPROBABILITY  Calculates any of the related (through Bayes rule) probabilities
%                 of a Gaussian Mixture Model (gmmDS) and a given dataset x.
%                 'prob_type' is a string indicating which of the four probability
%                 values are needed. These probabilities are:
%
%                      P(X|C) . P(C)                       likelihood . prior
%           P(C|X) = -----------------       posterior =  --------------------
%                          P(X)                                evidence
%
%            where C is the component classes (Gaussians) of the GMM and X is the data.
%
%   probability = gmmprobability(gmmDS, X, W)
%
%   INPUT
%          gmmDS         Gaussian mixture model data structure with the following fields
%            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
%            .dim        data dimension  [scalar]
%            .M          number of Gaussian component densities  [scalar]
%            .weights    mixing priors (component weights) [1-by-M vector]
%            .mu         M Gaussian component means (columns of matrix) [dim-by-M matrix]
%            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
%                        [dim-by-dim-by-M matrix]
%          X             buffer of N dim-by-1 data set vectors to be evaluated  [dim-by-N]
%          W             (optional) 1-by-N vector of sample weights. If specified, the sample
%                                   set will be weighted according to these weights.
%
%   OUTPUT
%          prior         The prior (without seeing data) probability of a component
%                        density generating any given data vector, i.e. P(C(i)).
%                        This is simply the same as the prior mixing weights,
%                        'gmmDS.weights'. [M-by-1 matrix]
%
%          likelihood    M-by-N martrix where the j,i-th entry is the likelihood
%                        of input column vector i (of X) conditioned on component
%                        density j, i.e. P(X(i)|C(j))
%
%          evidence      1-by-N matrix where the i-th entry is the total data probability
%                        for a given data vector X(i), i.e. P(X(i))=sum_over_all_j[P(X(i)|C(j))]
%
%          posterior     M-by-N matrix where the j,i-th entry is the posterior
%                        probability (after seeing the data) that a component
%                        density j has generated a specific data vector X(i), i.e.
%                        P(C(j)|X(i))   (class posterior probabilities)

    nOut = nargout;

    [dim, nov] = size(x);

    if (dim ~= gmmDataSet.dimension)
        error(' [ gmmProbability ] Data dimension and GMM model dimension is not the same.');
    end

    M        = gmmDataSet.M;                      % dumber of component densities
    mu       = gmmDataSet.mean;
    covar    = gmmDataSet.covariance;
    prior    = gmmDataSet.weights(:);             % prior mixing probabilities
    ones_nov = ones(nov,1);
    ones_M   = ones(M,1);

    if nOut > 1

        likelihood = zeros(M, nov);
        normfact = (2*pi) ^ (gmmDataSet.dimension / 2);

        switch gmmDataSet.covarianceType
            case {'full','diag'}            
                for k=1:M,
                    cmu = mu(:,k);
                    xx = x - cmu(:,ones_nov);
                    s = chol(covar(:,:,k))';
                    foo = s \ xx;
                    likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(s))));
                end

            case {'sqrt','sqrt-diag'}            
                for k=1:M,
                    cmu = mu(:,k);
                    xx = x - cmu(:,ones_nov);
                    s = covar(:,:,k);
                    foo = s \ xx;
                    likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(s))));
                end

            otherwise            
                error([' [ gmmProbability ] Unknown covariance type ', mix.cov_type]);

        end

    end

    likelihood = likelihood + 1e-99; % avoid zero

    % Calculate evidence
    if nOut > 2    
        if (nargin == 3)
            evidence = prior' * (likelihood ./ w(ones_M,:)); % weighted
        else
            evidence = prior'*likelihood; % non-weighted
        end

        evidence = evidence + 1e-99; % avoid zero
    end

    % Calculate posterior
    if nOut > 3    
        posterior = likelihood ./ ((1./prior)*evidence) + 1e-99;    
        posterior = posterior ./ rvecrep(sum(posterior,1),M);

    end
end