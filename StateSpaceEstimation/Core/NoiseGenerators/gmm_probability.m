function [ prior_mixing_probabilities, likelihood, evidence, posterior ] = gmm_probability(gmm_model, particles, weights )
    % gmm_probability. Calculates any of the related (through Bayes rule) probabilities of a Gaussian Mixture Model (gmmSet) and a given dataset (dataSet).
    %       Probabilities are:
    %                   P(X|C) . P(C)                   likelihood . prior
    %       P(C|X) = -----------------   posterior =  --------------------
    %                      P(X)                             evidence
    %
    %   where C is the component classes (Gaussians) of the GMM and X is the data.
    %
    %   [ prior, likelihood, evidence, posterior ] = gmm_probability( gmmSet, dataSet, evidenceWeights )
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
    output_arg_count = nargout;
    
    % number and dimension of input vectors
    [dimension, input_vect_count] = size(particles);
    
    if (dimension ~= gmm_model.dimension)
        error(' [ gmm_probability ] Data dimension and GMM model dimension is not the same.');
    end
    
    mixtureCount = gmm_model.mixtureCount;
    mu    = gmm_model.mean;
    covar = gmm_model.covariance;
    prior_mixing_probabilities = gmm_model.weights(:);
    
    ones_input_vectors  = ones(input_vect_count, 1);
    ones_mixtures       = ones(mixtureCount, 1);
    
    
    % Calculate likelihood
    if output_arg_count > 1
        likelihood = zeros(mixtureCount,input_vect_count);
        normfact = (2*pi)^(gmm_model.dimension / 2);
        
        for k = 1:mixtureCount
            switch gmm_model.covarianceType
                case {'full','diag'}
                    cov_x = chol(covar(:, :, k))';
                case {'sqrt','sqrt-diag'}
                    cov_x = covar(:, :, k);
                otherwise
                    error([' [ gmm_probability ] Unknown covariance type ', gmm_model.covarianceType]);
            end
            
            cmu = mu(:, k);
            x_centered = particles - cmu(:, ones_input_vectors);
            foo = cov_x \ x_centered;
            likelihood(k, :) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(cov_x))));
        end
    end
    likelihood = likelihood + 1e-99;
    
    
    % Calculate evidence
    if output_arg_count > 2
        if (nargin == 3)
            evidence = prior_mixing_probabilities' * (likelihood ./ weights(ones_mixtures, :));
        else
            evidence = prior_mixing_probabilities'*likelihood;
        end
        
        evidence = evidence + 1e-99;
    end
    
    
    % Calculate posterior
    if output_arg_count > 3
        posterior = likelihood ./ ((1./prior_mixing_probabilities)*evidence) + 1e-99;
        % normalize posterior
        posterior = posterior ./ row_vector_replicate( sum(posterior, 1), mixtureCount );
    end
end
