function noise = sampleGaussian(gausDataSet, n)
    % sampleGaussian Draw N samples from the Gaussian distribution (pdf) described by the
    %   Gaussian data structure 'gausDataSet'.
    %   noise = sampleGaussian(gausDataSet, n)
    %
    %   INPUT
    %          gausDataSet       Gaussian data structure with the following fields
    %             .covarianceType   (string)   covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'
    %             .dimension        (scalar)   dimension.
    %             .mean             (c-vector) mean vector  (dim-by-1)
    %             .covariance       (matrix)   covariance matrix of type covarianceType  (dimension-by-dimension)
    %          n                    (scalar)   number of samples to generate
    %   OUTPUT
    %       noise (matrix)   buffer of generated samples (dim-by-N)
    %
    %%
    switch gausDataSet.covarianceType
        case {'full','diag'}
            sqrtCov = chol(gausDataSet.covariance, 'lower');
        case {'sqrt','sqrt-diag'}
            sqrtCov = gausDataSet.covariance;
        otherwise
            error('[ sampleGaussian ] Unknown covariance type ');
    end
    
    noise = sqrtCov * randn(gausDataSet.dimension, n) + gausDataSet.mean(:, ones(n, 1));
end
