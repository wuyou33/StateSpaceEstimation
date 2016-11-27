function [ noise ] = sampleGaussian( gaussDataSet, count )
    % sampleGaussian. Draw N (count samples from the Gaussian distribution (pdf) described by the Gaussian data structure 'gaussDataSet'.
    %
    %   [ noise ] = sampleGaussian( gaussDataSet, count )
    %
    %   INPUT
    %          gaussDataSet          Gaussian data structure with the following fields
    %             .covarianceType        covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag', 'svd';
    %             .dimension             stochastic process dimension;
    %             .mean                  mean vector (dimension-by-1);
    %             .covariance            covariance matrix with type covarianceType  (dimension-by-dimension);
    %          count                 number of samples to generate.
    %
    %   OUTPUT
    %       noise   buffer of generated samples (dimension-by-count).
    %
    narginchk(2, 2);
    
    switch gaussDataSet.covarianceType
        case {'full', 'diag'}
            sqrtCov = chol(gaussDataSet.covariance, 'lower');
        case {'sqrt', 'sqrt-diag', 'svd'}
            sqrtCov = gaussDataSet.covariance;
        otherwise
            error('[ sampleGaussian::gaussDataSet ] Unknown covariance type ');
    end
    
    noise = sqrtCov * randn(gaussDataSet.dimension, count) + gaussDataSet.mean(:, ones(count, 1));
end
