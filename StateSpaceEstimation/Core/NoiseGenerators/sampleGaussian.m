function noise = sampleGaussian(gausDataSet, n)
%% GAUSSAMPLE  Draw N samples from the Gaussian distribution (pdf) described by the Gaussian data structure 'gausDataSet'.
%%
    noise = gausDataSet.covariance * randn(gausDataSet.dimension, n) + gausDataSet.mean(:, ones(n, 1));
end