function [ llh ] = likelihoodComboGaussian( noiseDataSet, noise, idxVec )
% Calculates the likelihood of sample 'noise', given the noise model NoiseDS. If the optional index
% vector 'idxVec' is specified, only those sub-noise sources are used. The 'noise' vector's dimension should
% concur with the implied total dimensionality of 'idxVec'

    if (nargin == 2)
        idxVec = 1:noiseDataSet.N;
    end

    numNS = length(idxVec);

    [~, nov] = size(noise);

    idxArr = noiseDataSet.idxArr;

    llh = ones(1,nov);

    for j=1:numNS,

        idx1 = idxArr(idxVec(j),1);
        idx2 = idxArr(idxVec(j),2);

        idxRange = idx1:idx2;

        dim = idx2-idx1+1;

        switch noiseDataSet.covarianceType

            case 'full'
                D  = det(noiseDataSet.covariance(idxRange,idxRange));
                iP = inv(noiseDataSet.covariance(idxRange,idxRange));
            case 'diag'
                D = prod(diag(noiseDataSet.covariance(idxRange,idxRange)));
                iP = diag(1./diag(noiseDataSet.covariance(idxRange,idxRange)));
            case 'sqrt'
                D = det(noiseDataSet.covariance(idxRange,idxRange))^2;
                iS = inv(noiseDataSet.covariance(idxRange,idxRange));
                iP = iS'*iS;
            case 'sqrt-diag'
                D = prod(diag(noiseDataSet.covariance(idxRange,idxRange)))^2;
                iP = diag(1./(diag(noiseDataSet.covariance(idxRange,idxRange)).^2));
            otherwise
                error(' [ likelihoodComboGaussian ] unknown cov_type.');
        end

        X = noise - cvecrep(noiseDataSet.mu(idxRange),nov);
        q = 1/sqrt((2*pi)^dim * D);

        llh = llh .* (q * exp(-0.5*diag(X'*iP*X)'));

    end

    llh = llh + 1e-99; % needed to avoid 0 likelihood (cause ill conditioning)
end
