function [gmm_set, log_evidence_buffer] = gmm_fit(samples, mixture_count, termination_threshold, covariance_type, check_cov, display, evidence_weights)
    % gmm_fit   Fit a Gaussian mixture model (GMM) with M components to dataset (gmmSet)
    % using an expectation maximization (EM) algorithm.
    %
    %   For more details see
    %       http://stats.stackexchange.com/questions/72774/numerical-example-to-understand-expectation-maximization
    %       http://stackoverflow.com/questions/11808074/what-is-an-intuitive-explanation-of-expectation-maximization-technique
    %       http://stackoverflow.com/questions/15513792/expectation-maximization-coin-toss-examples
    %
    %   [gmm_set, log_evidence_buffer] = gmm_fit(data, mixture_count, termination_threshold, covariance_type, check_cov, display, evidence_weights)
    %
    %   INPUT
    %       samples                 Dataset of N samples (column vectors);
    %       mixture_count           Number of Gaussian mixture component densities OR a pre-initialized GMM data structure;
    %       termination_threshold 	Termination threshold 0 < tt < 1 (if change in log likelihood falls below this value, the EM algorithm terminates.)
    %                                   OR if tt is a [1-by-2 vector] the first component is the termination threshold and the second component is the maximum
    %                                   number of iterations allowed for the EM, i.e. thresholdtt iterations];
    %       covariance_type       	Covariance type 'full', 'diag', 'sqrt', 'sqrt-diag';
    %       check_cov       	    <optional> Covariance check flag: If this flag is true, then covariance matrix is reset to its original value when any of its
    %                                   singular values are too small;
    %       display                 log temporary results;
    %       evidence_weights      	<optional> vector of sample weights used to do a weighted EM fit to the data.
    %
    %   OUTPUT
    %       gmm_set              - Gaussian mixture model data structure with the following fields;
    %           .covarianceType - Covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag';
    %           .dimension      - Data dimension;
    %           .mixtureCount   - Number of Gaussian component densities;
    %           .weights        - Mixing priors (component weights);
    %           .mean           - Gaussian component means;
    %           .covariance     - Covariance matrices of Gaussian components (must comply with cov_type);
    %       log_evidence_buffer - log evidence buffer (sum of log evidence of data as a function of EM iteration).
    %
    switch nargin
        case 2
            termination_threshold = [0.1 100];
            covariance_type = 'full';
            check_cov = 0;
            display = 0;
        case 3
            covariance_type = 'full';
            check_cov = 0;
            display = 0;
        case 4
            check_cov = 0;
            display = 0;
        case 5
            display = 0;
        case {6,7}
        otherwise
            error(' [ gmm_fit ] Incorrect number of input arguments.');
    end
    
    [dim, data_size] = size(samples);
    
    % Sort out the options
    if (length(termination_threshold) == 2)
        iterations_count = termination_threshold(2);
    else
        iterations_count = 100;
    end
    
    store = 0;
    if (nargout > 1)
        store = 1;
        log_evidence_buffer = zeros(1, iterations_count);
    end
    
    test = 0; % Test log likelihood for termination
    if (termination_threshold(1) > 0.0)
        test = 1;
        log_evidence_buffer = zeros(1, iterations_count);
    end
    
    if (isnumeric(mixture_count))
        gmm_set.covarianceType = covariance_type;
        gmm_set.dimension = dim;
        gmm_set.mixtureCount = mixture_count;
        gmm_set.weights = ones(dim, mixture_count);
        gmm_set.mean = zeros(dim, mixture_count);
        gmm_set.covariance = repmat(eye(dim), [1 1 mixture_count]);
        gmm_set = gmm_initialize(gmm_set, samples, 5);
    else
        gmm_set = mixture_count;
        mixture_count = gmm_set.mixtureCount;
    end
    
    % Ensure that covariances don't collapse
    if check_cov
        min_cov = eps;
        min_cov_sqrt = sqrt(eps);
        init_covars = gmm_set.covariance;
    end
    
    eold = -Inf;
    ones_dim = ones(1, dim);
    ones_data_size = ones(data_size, 1);
    
    for n = 1 : iterations_count
        % Calculate posteriors based on old parameters
        if (nargin == 7)
            [~, ~, evidence, posterior] = gmm_probability(gmm_set, samples, evidence_weights);
        else
            [~, ~, evidence, posterior] = gmm_probability(gmm_set, samples);
        end
        
        % Calculate error value if needed
        if (display || store || test)
            % Error value is negative log likelihood of data evidence
            err_val = sum(log(evidence));
            if store
                log_evidence_buffer(n) = err_val;
            end
            
            if display
                fprintf(1, 'Cycle %4d  Evidence %11.6f\n', n, err_val);
            end
            
            if test
                if (n > 1 && abs((err_val - eold)/eold) < termination_threshold(1))
                    log_evidence_buffer=log_evidence_buffer(1:n);
                    return;
                else
                    eold = err_val;
                end
            end
        end
        
        
        % Adjust the new estimates for the parameters
        new_pr = (sum(posterior, 2))';
        new_c =  (posterior * samples')';
        
        % Now move new estimates to old parameter vectors
        gmm_set.weights = new_pr / data_size;
        new_mean = new_c ./ new_pr(ones_dim, :);
        gmm_set.mean = new_mean;
        
        for j = 1:mixture_count
            new_mean_j = new_mean(:, j);
            diffs = samples - new_mean_j(:, ones_data_size);
            
            switch covariance_type
                case 'full'
                    t_post = sqrt(posterior(j, :));
                    diffs = diffs .* t_post(ones_dim, :);
                    gmm_set.covariance(:, :, j) = (diffs*diffs') / new_pr(j);
                    
                    if check_cov
                        if min(svd(gmm_set.covariance(:, :, j))) < min_cov
                            gmm_set.covariance(:, :, j) = init_covars(:, :, j);
                        end
                    end
                case {'sqrt','sqrt-diag'}
                    t_post = ( 1 / sqrt(new_pr(j)) ) * sqrt(posterior(j, :));
                    diffs = diffs .* t_post(ones_dim, :);
                    [~, tcov] = qr(diffs', 0);
                    gmm_set.covariance(:, :, j) = tcov';
                    
                    if check_cov
                        if min(abs(diag(gmm_set.covariance(:, :, j)))) < min_cov_sqrt
                            gmm_set.covariance(:, :, j) = init_covars(:, :, j);
                        end
                    end
                case 'diag'
                    t_post = posterior(j, :);
                    dd = sum( (diffs.*diffs) .* t_post(ones_dim,:), 2 ) /new_pr(j);
                    
                    if check_cov && min(dd) < min_cov
                        gmm_set.covariance(:, :, j) = init_covars(:, :, j);
                    else
                        gmm_set.covariance(:, :, j) = diag(dd);
                    end
                otherwise
                    error(['Unknown covariance type ', covariance_type]);
            end
        end
    end
end
