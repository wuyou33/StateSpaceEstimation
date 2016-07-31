function [ points, edges ] = gmmSample( gmmSet, count )
    % gmmSample  Draw N samples from the Gaussian mixture model (GMM) described by the GMM data structure 'gmmSet'.
    %
    %   [x, index] = gmmSample(gmmSet, n)
    %
    %   INPUT
    %       gmmDS               - Gaussian mixture model data structure with the following fields;
    %           .covarianceType - covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag';
    %           .dimension      - data dimension;
    %           .mixtureCount   - number of Gaussian component densities;
    %           .weights        - mixing priors (component weights);
    %           .mean           - Gaussian component means (columns of matrix);
    %           .covariance     - covariance matrices of Gaussian components;
    %       count               - number of samples to generate.
    %   OUTPUT
    %       points - buffer of N samples drawn from the GMM;
    %       edges  - component index of samples.
    %
    %%
    narginchk(2, 2);
    xrnd = rand(1, count);
    
    % draw component indices according to prior
    [n, edges, bin] = histcounts(xrnd, cumsum([0; gmmSet.weights(:)]));
    points = zeros(gmmSet.dimension, count);
    
    switch gmmSet.covarianceType
        case {'full','diag'}
            for i = 1:gmmSet.mixtureCount
                idx = (find(bin == i));
                points(:, idx) = chol(gmmSet.covariance(:, :, i), 'lower') * randn(gmmSet.dimension, n(i));
            end
        case {'sqrt','sqrt-diag'}
            for i = 1:gmmSet.mixtureCount
                idx = (find(bin == i));
                points(:, idx) = gmmSet.covariance(:, :, i) * randn(gmmSet.dimension, n(i));
            end
        otherwise
            error(' [ gmmSample ] Unknown covariance type');
    end
    
    points = points + gmmSet.mean(:, bin);
end
