function [ prior, likelihood, evidence, posterior ] = gmmProbability( gmmDS, X, W )
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
    %           .covarianceType    covariance matrix type ('full' , 'diag' , 'sqrt' , 'sqrt-diag');
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
    Nout = nargout;
    
    [dim,nov] = size(X);              % number and dimension of input vectors
    
    if (dim~=gmmDS.dimension)
        error(' [ gmmprobability ] Data dimension and GMM model dimension is not the same.');
    end
    
    M = gmmDS.mixtureCount;                      % dumber of component densities
    mu    = gmmDS.mean;                 % component means
    covar = gmmDS.covariance;                % component covariance matrices
    
    prior = gmmDS.weights(:);        % prior mixing probabilities
    
    ones_nov = ones(nov,1);
    ones_M   = ones(M,1);
    
    %--- Calculate likelihood
    if Nout > 1
        
        likelihood = zeros(M,nov);        % preallocate component likelihood matrix
        normfact = (2*pi)^(gmmDS.dimension/2);  % component density normalizing factor
        
        switch gmmDS.covarianceType             % calculate per component likelihood
            
            case {'full','diag'}
                
                for k=1:M,
                    cmu = mu(:,k);
                    XX = X - cmu(:,ones_nov);
                    S = chol(covar(:,:,k))';
                    foo = S \ XX;
                    likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(S))));
                end
                
            case {'sqrt','sqrt-diag'}
                
                for k=1:M,
                    cmu = mu(:,k);
                    XX = X - cmu(:,ones_nov);
                    S = covar(:,:,k);
                    foo = S \ XX;
                    likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(S))));
                end
                
            otherwise
                
                error([' [ gmmprobability ] Unknown covariance type ', mix.covarianceType]);
                
        end
        
    end
    
    likelihood = likelihood + 1e-99;
    
    
    %--- Calculate evidence
    if Nout > 2
        
        if (nargin == 3)
            evidence = prior' * (likelihood ./ W(ones_M,:));  % weighted
        else
            evidence = prior'*likelihood;                     % non-weighted
        end
        
        evidence = evidence + 1e-99;
        
    end
    
    
    %--- Calculate posterior
    if Nout > 3
        
        posterior = likelihood ./ ((1./prior)*evidence) + 1e-99;
        % normalize
        posterior = posterior ./ rvecrep(sum(posterior,1),M);
        
    end
end
