function [gmmDS, leb] = gmm_fit(data, mixtureCnt, terminationThreshold, covarianceType, check_cov, evidenceWeights)
    % gmm_fit   Fit a Gaussian mixture model (GMM) with M components to dataset (gmmSet)
    % using an expectation maximization (EM) algorithm.
    %
    %   For more details see
    %       http://stats.stackexchange.com/questions/72774/numerical-example-to-understand-expectation-maximization
    %       http://stackoverflow.com/questions/11808074/what-is-an-intuitive-explanation-of-expectation-maximization-technique
    %       http://stackoverflow.com/questions/15513792/expectation-maximization-coin-toss-examples
    %
    %   [ gmmSet ] = gaussMixtureModelFit( data, mixtureCnt, terminationThreshold, covarianceType, evidenceWeights )
    %
    %   INPUT
    %       data                 	Dataset of N samples (column vectors);
    %       mixture              	Number of Gaussian mixture component densities OR a pre-initialized GMM data structure;
    %       terminationThreshold 	Termination threshold 0 < tt < 1 (if change in log likelihood falls below this value, the EM algorithm terminates.)
    %                                   OR if tt is a [1-by-2 vector] the first component is the termination threshold and the second component is the maximum
    %                                   number of iterations allowed for the EM, i.e. thresholdtt iterations];
    %       covarianceType       	Covariance type 'full', 'diag', 'sqrt', 'sqrt-diag';
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
    %
    %   INPUT
    %          X            Dataset of N samples (column vectors) [dim-by-N matrix]
    %          M            Number of Gaussian mixture component densities [scalar]
    %                       OR a pre-initialized GMM data structure (gmmDS)
    %          tt           Termination threshold 0 < tt < 1 (if % change in log likelihood
    %                       falls below this value, the EM algorithm terminates.) [scalar]
    %                       OR if tt is a [1-by-2 vector] the first component is the termination
    %                       threshold and the second component is the maximum number of iterations
    %                       allowed for the EM, i.e. tt = [tt max_iterations]
    %          cov_type     Covariance type 'full','diag','sqrt','sqrt-diag' [string]
    %          check_cov    (optional) Covariance check flag : If this flag is set, a covariance
    %                       matrix is reset to its original value when any of its singular values
    %                       are too small (less than MIN_COVAR which has the value eps). With the
    %                       default value of 0 no action is taken.
    %          display      (optional) Display results of training if this is set to 1 (default=0)
    %          W            (optional) 1xN vector of sample weights used to do a weighted EM fit to
    %                                  the data.
    %
    %   OUTPUT
    %          gmmDS         Gaussian mixture model data structure with the following fields
    %            .cov_type   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'    [string]
    %            .dim        data dimension  [scalar]
    %            .M          number of Gaussian component densities  [scalar]
    %            .weights    mixing priors (component weights) [1-by-M matrix]
    %            .mu         N Gaussian component means (columns of matrix) [dim-by-N matrix]
    %            .cov        covariance matrices of Gaussian components (must comply with .cov_type)
    %                        [dim-by-dim-by-N matrix]
    %          leb           log evidence buffer (sum of log evidence of data as a function of EM iteration)
    
    narginchk(2, 6);
    switch nargin
        case 2
            terminationThreshold = [0.1 100];
            covarianceType = 'full';
            check_cov = 0;
        case 3
            covarianceType = 'full';
            check_cov = 0;
        case 4
            check_cov = 0;
    end
    
    display = 0;
    [dim, Ndata] = size(data);             % get dimensions of dataset
    
    % Sort out the options
    if (length(terminationThreshold)==2)                  % number of EM iterations
        niters = terminationThreshold(2);
    else
        niters = 100;
    end
    
    
    store = 0;
    if (nargout > 1)
        store = 1;    % Store the evidence values to return them
        leb = zeros(1, niters);
    end
    
    test = 0;
    if terminationThreshold(1) > 0.0
        test = 1;                         % Test log likelihood for termination
        leb = zeros(1, niters);
    end
    
    
    if isnumeric(mixtureCnt)                    % Initialize GMM from data if needed
        gmmDS.covarianceType = covarianceType;
        gmmDS.dimension = dim;
        gmmDS.mixtureCount = mixtureCnt;
        gmmDS.weights = ones(dim,mixtureCnt);
        gmmDS.mean = zeros(dim,mixtureCnt);
        gmmDS.covariance = repmat(eye(dim),[1 1 mixtureCnt]);   % initialize GMM centroids and their covariances
        gmmDS = gmm_initialize(gmmDS, data, 5);     % using Netlab's KMEANS algorithm.
    else
        gmmDS = mixtureCnt;
        mixtureCnt = gmmDS.mixtureCount;
    end
    
    
    if check_cov                        % Ensure that covariances don't collapse
        MIN_COVAR = eps;                  % Minimum singular value of covariance matrix
        MIN_COVAR_SQRT = sqrt(eps);
        init_covars = gmmDS.covariance;
    end
    
    
    eold = -Inf;
    
    ones_dim = ones(1,dim);
    ones_Ndata = ones(Ndata,1);
    
    % Main loop of algorithm
    for n = 1:niters
        
        % Calculate posteriors based on old parameters
        if (nargin == 7)
            [~, ~, evidence, posterior] = gmmProbability(gmmDS, data, evidenceWeights);
        else
            [~, ~, evidence, posterior] = gmmProbability(gmmDS, data);
        end
        
        % Calculate error value if needed
        if (display || store || test)
            % Error value is negative log likelihood of data evidence
            e = sum(log(evidence));
            if store
                leb(n) = e;
            end
            if display
                fprintf(1, 'Cycle %4d  Evidence %11.6f\n', n, e);
            end
            if test
                if (n > 1 && abs((e - eold)/eold) < terminationThreshold(1))
                    leb=leb(1:n);
                    return;
                else
                    eold = e;
                end
            end
        end
        
        
        % Adjust the new estimates for the parameters
        new_pr = (sum(posterior, 2))';
        new_c =  (posterior * data')';
        
        % Now move new estimates to old parameter vectors
        gmmDS.weights = new_pr / Ndata;
        new_mu = new_c ./ new_pr(ones_dim,:);
        gmmDS.mean = new_mu;
        
        switch covarianceType
            
            case 'full'
                for j = 1:mixtureCnt
                    tmu = new_mu(:,j);
                    diffs = data - tmu(:,ones_Ndata);
                    tpost = sqrt(posterior(j,:));
                    diffs = diffs .* tpost(ones_dim,:);
                    gmmDS.covariance(:,:,j) = (diffs*diffs')/new_pr(j);
                end
                if check_cov
                    % Ensure that no covariance is too small
                    for j = 1:mixtureCnt
                        if min(svd(gmmDS.covariance(:,:,j))) < MIN_COVAR
                            gmmDS.covariance(:,:,j) = init_covars(:,:,j);
                        end
                    end
                end
                
            case {'sqrt','sqrt-diag'}
                for j = 1:mixtureCnt
                    tmu = new_mu(:,j);
                    diffs = data - tmu(:,ones_Ndata);
                    tpost = (1/sqrt(new_pr(j))) * sqrt(posterior(j,:));
                    diffs = diffs .* tpost(ones_dim,:);
                    [~,tcov] = qr(diffs',0);
                    gmmDS.covariance(:,:,j) = tcov';
                end
                if check_cov
                    % Ensure that no covariance is too small
                    for j = 1:mixtureCnt
                        if min(abs(diag(gmmDS.covariance(:,:,j)))) < MIN_COVAR_SQRT
                            gmmDS.covariance(:,:,j) = init_covars(:,:,j);
                        end
                    end
                end
                
                
            case 'diag'
                for j = 1:mixtureCnt
                    tmu = new_mu(:,j);
                    diffs = data - tmu(:,ones_Ndata);
                    tpost = posterior(j,:);
                    dd = sum((diffs.*diffs).*tpost(ones_dim,:), 2)/new_pr(j);
                    if check_cov
                        if min(dd) < MIN_COVAR
                            gmmDS.covariance(:,:,j) = init_covars(:,:,j);
                        else
                            gmmDS.covariance(:,:,j) = diag(dd);
                        end
                    else
                        gmmDS.covariance(:,:,j) = diag(dd);
                    end
                end
                
                
            otherwise
                error(['Unknown covariance type ', covarianceType]);
        end
        
    end
end
