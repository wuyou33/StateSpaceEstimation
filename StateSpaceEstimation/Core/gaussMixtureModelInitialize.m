function [ gmmSet ] = gaussMixtureModelInitialize( gmmSet, dataSet, iterationNumber, arbitraryWidth )
    % gaussMixtureModelInitialize Initialises Gaussian mixture model (GMM) from data
    %
    %   [ gmmSet ] = gaussMixtureModelInitialize( gmmSet, dataSet, iterationNumber, arbitraryWidth )
    %
    %   INPUT
    %       gmmSet              - Gaussian mixture model data structure with the following fields;
    %           .covarianceType - covariance matrix type ('full' , 'diag' , 'sqrt' , 'sqrt-diag');
    %           .dimension      - data dimension;
    %           .mixtureCount   - number of Gaussian component densities;
    %           .weights        - component weights;
    %           .mean           - Gaussian component means;
    %           .covariance     - covariance matrices of Gaussian components (must comply with covarianceType);
    %       dataSet             - dataset of mixtureCount samples;
    %       iterationNumber     - <optional> maximum number of iterations (default value is 200).
    %       arbitraryWidth      - <optional>  Arbitrary width used if variance collapses to zero: make it 'large' 
    %                                         so that centre is responsible for a reasonable number of points.
    %
    %   OUTPUT
    %       gmmSet  - initilized / updated GMM data structure.
    %%
    narginchk(2, 4);
    %%
    if nargin < 3; iterationNumber = 200; end
    if nargin < 4; arbitraryWidth = 1.0; end
    
    stateDim = size(dataSet, 1);
    dataSet = dataSet';
    
    % Use kmeans algorithm to initialise the centroids from the data
    options.logError = 0;
    options.initializeFromData = 1;
    options.absPrecision = 1e-5;
    options.precision = 1e-5;
    options.iterationNumber = iterationNumber;
    [mean, ~, post] = kmeans2(gmmSet.mean', dataSet, options);
    gmmSet.mean = mean';    
    
    clusterSizes = max(sum(post, 1), 1);
    gmmSet.weights = clusterSizes / sum(clusterSizes);
    
    fixCov = arbitraryWidth*eye(stateDim);    
    for j = 1 : gmmSet.mixtureCount
        idx = (find(post(:, j)));
        centroidPoints = dataSet(idx, :);
        sizec = size(centroidPoints, 1);
        mu = mean(j, :);
        diffs = centroidPoints - rvecrep(mu, sizec);
        
        switch gmmSet.covarianceType
            case 'full'
                gmmSet.covariance(:, :, j) = (diffs'*diffs) / sizec;
                
                if rank(gmmSet.covariance(:, :, j)) < stateDim
                    gmmSet.covariance(:, :, j) = gmmSet.covariance(:, :, j) + fixCov;
                end
            case 'diag'
                cov = sum((diffs.*diffs), 1)/sizec;
                cov = cov + arbitraryWidth*(cov < eps);
                gmmSet.covariance(:, :, j) = diag(cov);
            case 'sqrt'
                cov = (diffs'*diffs)/sizec;
                if rank(cov) < gmmSet.dimension;
                    cov = cov + fixCov;
                end
                gmmSet.covariance(:, :, j) = chol(cov, 'lower');
            case 'sqrt-diag'
                cov = sum((diffs.*diffs), 1) / sizec;
                cov = cov + arbitraryWidth*(cov < eps);
                gmmSet.covariance(:,:,j) = diag(sqrt(cov));
            otherwise
                error([' [ gaussMixtureModelInitialize ] Unknown covariance type ', gmmSet.covarianceType]);
        end
    end
end
