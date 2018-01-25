function gmm_set = gmm_initialize(gmm_set, particles, max_iterations)
    %   gmm_initialize. Initialises Gaussian mixture model (GMM) from input dataSet with specified parameters.
    %
    %   gmm_set = gmm_initialize(gmm_set, particles, max_iterations)
    %
    %   INPUT
    %       gmm_set              	Gaussian mixture model data structure with the following fields;
    %           .covarianceType         covariance matrix type ('full' , 'diag' , 'sqrt' , 'sqrt-diag');
    %           .dimension              data dimension;
    %           .mixtureCount           number of Gaussian component densities;
    %           .weights                component weights;
    %           .mean                   Gaussian component means;
    %           .covariance             covariance matrices of Gaussian components (must comply with covarianceType);
    %       particles             	dataset of mixture count samples;
    %       max_iterations          <optional> maximum number of iterations (default = 100).
    %
    %   OUTPUT
    %       gmm_set     initilized / updated GMM data structure.
    %
    if nargin < 3
        max_iterations = 100;
    end
    
    [state_dim, ~] = size(particles);    
    mixture_count   = gmm_set.mixtureCount;
    particles = particles';
    
    % Arbitrary width used if variance collapses to zero: make it 'large' so
    % that centre is responsible for a reasonable number of points.
    arbitrary_width_gmm = 1.0;
    
    % Use kmeans algorithm to initialise the centroids from the data
    options = foptions;
    options(1) = -1; % don't display warnings
    options(14) = max_iterations; % Just use max_iterations iterations of k-means in initialisation
    options(5)  = 1; % initilize centroids and their covariances from data
    [mu, ~, post] = kmeans3(gmm_set.mean', particles, options);
    
    gmm_set.mean = mu';    
    % set priors depending on number of points in each cluster
    cluster_sizes = max(sum(post, 1), 1); % Make sure that no prior is zero
    gmm_set.weights = cluster_sizes / sum(cluster_sizes); % normalise priors
    
    fix_cov = arbitrary_width_gmm * eye(state_dim);
    for j = 1:mixture_count
        % Pick out data points belonging to this centre
        centres = particles(find(post(:, j)), :);        
        size_centres = size(centres, 1);
        mean_center = mu(j, :);
        diffs = centres - mean_center(ones(1, size_centres), :);
        
        switch gmm_set.covarianceType
            case 'full'
                gmm_set.covariance(:, :, j) = (diffs' * diffs) / size_centres;
                % Add arbitrary_width_gmm * Identity to rank-deficient covariance matrices
                if rank(gmm_set.covariance(:, :, j)) < state_dim
                    gmm_set.covariance(:, :, j) = gmm_set.covariance(:,:,j) + fix_cov;
                end
            case 'diag'
                diag_cov_mix = sum(diffs.*diffs, 1) / size_centres;
                % Replace small entries by arbitrary_width_gmm value
                diag_cov_mix = diag_cov_mix + arbitrary_width_gmm * (diag_cov_mix < eps);
                gmm_set.covariance(:, :, j) = diag(diag_cov_mix);
            case 'sqrt'
                sqrt_cov_mix = (diffs' * diffs) / size_centres;
                % Add arbitrary_width_gmm*Identity to rank-deficient covariance matrices
                if rank(sqrt_cov_mix) < gmm_set.dimension
                    sqrt_cov_mix = sqrt_cov_mix + fix_cov;
                end
                gmm_set.covariance(:, :, j) = chol(sqrt_cov_mix, 'lower');
            case 'sqrt-diag'
                diag_cov_mix = sum(diffs.*diffs, 1) / size_centres;
                % Replace small entries by arbitrary_width_gmm value
                diag_cov_mix = diag_cov_mix + arbitrary_width_gmm * (diag_cov_mix < eps);
                gmm_set.covariance(:, :, j) = diag(sqrt(diag_cov_mix));
            otherwise
                error([' [ gmm_initialize ] Unknown covariance type ', gmm_set.covarianceType]);
        end
    end
end
