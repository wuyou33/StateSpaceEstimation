function [ x, comp ] = gmmSample( gmmDataSet, n )
% GMMSAMPLE  Draw N samples from the Gaussian mixture model (GMM) described by the GMM data structure 'gmmDataSet'.
%%

    dim   = gmmDataSet.dim;                         % random vector dimension
    Ncomp = gmmDataSet.M;                           % number of component densities
    w     = gmmDataSet.weights(:);                  % prior mixing probabilities
    mu    = gmmDataSet.mu;                          % component means
    cov   = gmmDataSet.cov;                         % component covariance matrices
    u     = rand(1,n);

    [nc, comp] = histcounts(u, cumsum([0; w]));      % draw component indices according to prior probabilities specified in gmmDS.weights
                                                     % nc = number of samples in each component comp = index vector

    x = zeros(dim, n);

    switch gmmDataSet.cov_type
        case {'full','diag'}

            for k = 1:Ncomp,
                idx = find(comp == k);
                x(:, idx) = chol(cov(:,:, k))' * randn(dim, nc(k));
            end

            %----------------------------------------------------------------------
        case {'sqrt','sqrt-diag'}

            for k= 1:Ncomp,
                idx = find(comp == k);
                x(:, idx) = cov(:,:, k) * randn(dim, nc(k));
            end

        otherwise
            error(' [ gmmSample ] Unknown covariance type!');

    end

    x = x + mu(:,comp);     % add in means
end