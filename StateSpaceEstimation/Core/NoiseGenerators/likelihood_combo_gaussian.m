function [ llh ] = likelihood_combo_gaussian(noiseDataSet, noise, idxVec)
    % likelihood_combo_gaussian. Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    % If the optional index vector 'idxVec' is specified, only those sub-noise sources are used.
    % The 'noise' vector's dimension should concur with the implied total dimensionality of 'idxVec'
    %
    %   [ llh ] = likelihood_combo_gaussian(noiseDataSet, noise, idxVec)
    %
    %   INPUT
    %       noiseDataSet    structure wich fully describe stochastic process;
    %       noise           sample vector of stochastic process (concreate vector of stochastic process);
    %       idxVec          <<optional>> specifies what sources should be used to calculate likelihood.
    %
    %   OUTPUT
    %       llh    likelihood of sample 'noise'.
    %
    narginchk(2, 3);
    
    if nargin == 2
        idxVec = 1 : noiseDataSet.N;
    end
    
    count = length(idxVec);
    nov = size(noise, 2);
    
    idxArr = noiseDataSet.idxArr;
    
    llh = ones(1, nov);
    
    for j = 1:count
        idx1 = idxArr(idxVec(j),1);
        idx2 = idxArr(idxVec(j),2);
        
        idxRange = idx1:idx2;
        
        dim = idx2-idx1+1;
        
        switch noiseDataSet.covarianceType
            case 'full'
                determinant  = det(noiseDataSet.covariance(idxRange, idxRange));
                invertCov = inv(noiseDataSet.covariance(idxRange, idxRange));
            case 'diag'
                determinant = prod(diag(noiseDataSet.covariance(idxRange, idxRange)));
                invertCov = diag(1./diag(noiseDataSet.covariance(idxRange, idxRange)));
            case 'sqrt'
                determinant = det(noiseDataSet.covariance(idxRange, idxRange))^2;
                invertSqrt = inv(noiseDataSet.covariance(idxRange, idxRange));
                invertCov = invertSqrt'*invertSqrt;
            case 'sqrt-diag'
                determinant = prod(diag(noiseDataSet.covariance(idxRange, idxRange)))^2;
                invertCov = diag(1./(diag(noiseDataSet.covariance(idxRange, idxRange)).^2));
            otherwise
                error('[ likelihood_combo_gaussian ] unknown covarianceType.');
        end
        
        X = noise - column_vector_replicate(noiseDataSet.mean(idxRange), nov);
        q = 1/sqrt((2*pi)^dim * determinant);
        
        llh = llh .* (q * exp(-0.5*diag(X'*invertCov*X)'));
        
    end
    
    % needed to avoid 0 likelihood
    llh = llh + 1e-99;
end
