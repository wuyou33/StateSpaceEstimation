function gmmSet = gmm_initialize(gmmSet, dataSet, maxIterations)
    %   gmm_initialize. Initialises Gaussian mixture model (GMM) from input dataSet with specified parameters.
    %
    %   [ gmmSet ] = gmm_initialize( gmmSet, dataSet, iterationNumber, arbitraryWidth )
    %
    %   INPUT
    %       gmmSet              	Gaussian mixture model data structure with the following fields;
    %           .covarianceType         covariance matrix type ('full' , 'diag' , 'sqrt' , 'sqrt-diag');
    %           .dimension              data dimension;
    %           .mixtureCount           number of Gaussian component densities;
    %           .weights                component weights;
    %           .mean                   Gaussian component means;
    %           .covariance             covariance matrices of Gaussian components (must comply with covarianceType);
    %       dataSet             	dataset of mixtureCount samples;
    %       maxIterations     	<optional> maximum number of iterations (default = 100).
    %
    %   OUTPUT
    %       gmmSet  initilized / updated GMM data structure.
    %
    
    if nargin < 3; maxIterations = 100; end
    
    [xdim, ~] = size(dataSet);
    dataSet=dataSet';
    
    % Arbitrary width used if variance collapses to zero: make it 'large' so
    % that centre is responsible for a reasonable number of points.
    GMM_WIDTH = 1.0;
    
    % Use kmeans algorithm to initialise the centroids from the data
    options = foptions;
    options(1) = -1;       % don't display warnings
    options(14) = maxIterations;    % Just use 5 iterations of k-means in initialisation
    options(5)  = 1;       % initilize centroids and their covariances from data
    [mu, ~, post] = kmeans3(gmmSet.mean', dataSet, options);       % call Netlab k-means algorithm
    
    gmmSet.mean = mu';   % convert from Netlab format to ReBEL format
    
    % Set priors depending on number of points in each cluster
    cluster_sizes = max(sum(post, 1), 1);             % Make sure that no prior is zero
    gmmSet.weights = cluster_sizes/sum(cluster_sizes); % Normalise priors
    
    
    fixCov = GMM_WIDTH*eye(xdim);
    
    switch gmmSet.covarianceType
        case 'full'
            for j = 1:gmmSet.mixtureCount
                % Pick out data points belonging to this centre
                c = dataSet(~(post(:, j)),:);
                sizec = size(c,1);
                tmu = mu(j,:);
                diffs = c - tmu(ones(1,sizec),:);
                gmmSet.cov(:,:,j) = (diffs'*diffs)/sizec;
                % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
                if rank(gmmSet.covariance(:,:,j)) < xdim
                    gmmSet.covariance(:,:,j) = gmmSet.covariance(:,:,j) + fixCov;
                end
            end
            
        case 'diag'
            for j = 1:gmmSet.mixtureCount
                % Pick out data points belonging to this centre
                c = dataSet(find(post(:, j)),:);
                sizec = size(c,1);
                tmu = mu(j,:);
                diffs = c - tmu(ones(1,sizec),:);
                d = sum((diffs.*diffs), 1)/sizec;
                % Replace small entries by GMM_WIDTH value
                d = d + GMM_WIDTH*(d<eps);
                gmmSet.covariance(:,:,j) = diag(d);
            end
            
        case 'sqrt'
            for j = 1:gmmSet.mixtureCount
                % Pick out data points belonging to this centre
                c = dataSet(find(post(:, j)),:);
                sizec = size(c,1);
                tmu = mu(j,:);
                diffs = c - tmu(ones(1,sizec),:);
                cov = (diffs'*diffs)/sizec;
                % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
                if rank(cov) < gmmSet.dimension
                    cov = cov + fixCov;
                end
                gmmSet.covariance(:,:,j) = chol(cov)';
            end
            
        case 'sqrt-diag'
            for j = 1:gmmSet.mixtureCount
                % Pick out data points belonging to this centre
                c = x(find(post(:, j)),:);
                sizec = size(c,1);
                tmu = mu(j,:);
                diffs = c - tmu(ones(1,sizec),:);
                d = sum((diffs.*diffs), 1)/sizec;
                % Replace small entries by GMM_WIDTH value
                d = d + GMM_WIDTH*(d<eps);
                gmmSet.covariance(:,:,j) = diag(sqrt(d));
            end
            
        otherwise
            error([' [ gmm_initialize ] Unknown covariance type ', gmmSet.covarianceType]);
    end
end
