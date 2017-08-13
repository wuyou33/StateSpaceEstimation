function [ gmmSet ] = gaussMixtureModelInitialize( gmmSet, dataSet, iterationNumber, arbitraryWidth, visualize )
    % gaussMixtureModelInitialize. Initialises Gaussian mixture model (GMM) from input dataSet with specified parameters.
    %
    %   [ gmmSet ] = gaussMixtureModelInitialize( gmmSet, dataSet, iterationNumber, arbitraryWidth )
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
    %       iterationNumber     	<optional> maximum number of iterations (default value is 200);
    %       arbitraryWidth      	<optional> arbitrary width used if variance collapses to zero: make it 'large'
    %                                  so that centre is responsible for a reasonable number of points;
    %       visualize               <optional> visualize clusters.
    %
    %   OUTPUT
    %       gmmSet  initilized / updated GMM data structure.
    %%
    narginchk(2, 5);
    
    if (nargin < 3)
        iterationNumber = 200;
    end
    
    if (nargin < 4)
        arbitraryWidth = 1.0;
    end
    
    if (nargin < 5)
        visualize = 0;
    end
    
    stateDim = size(dataSet, 1);
    dataSet  = dataSet';
    
    % Use kmeans algorithm to initialise the centroids from the data
    options.logError = 0;
    options.initializeFromData = 1;
    options.absPrecision = 1e-5;
    options.precision = 1e-5;
    options.iterationNumber = iterationNumber;
    options.visualize = visualize;
    
    [mean, ~, post] = kmeans2(gmmSet.mean', dataSet, options);
    gmmSet.mean = mean';
    
    clusterCnt = max(sum(post, 1), 1);
    gmmSet.weights = clusterCnt / sum(clusterCnt);
    fixCov = arbitraryWidth * eye(stateDim);
    
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
                error(['[ gaussMixtureModelInitialize::gmmSet ] Unknown covariance type ', gmmSet.covarianceType]);
        end
    end
    
    if (visualize)
        warning(' [ gaussMixtureModelInitialize ] visualizatoin is not correct!! Do not use...');
        startIdx = 1;
        endIdx   = 3;
        
        while (startIdx < stateDim)
            plotClaster(dataSet(:, startIdx : endIdx), length(clusterCnt));
            
            startIdx = startIdx + 3;
            
            if (endIdx + 3 < stateDim)
                endIdx = endIdx + 3;
            else
                endIdx = stateDim;
            end            
        end
    end
end
