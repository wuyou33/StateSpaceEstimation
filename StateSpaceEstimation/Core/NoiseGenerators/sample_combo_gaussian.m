function [ noise ] = sample_combo_gaussian( noiseDataSet, count )
    % sample_combo_gaussian. Generate N (count) samples of a noise source specified by the noiseDataSet data structure (mixture of Gaussian stochastic processes).
    %
    %   [ noise ] = sample_combo_gaussian( noiseDataSet, count )
    %
    %   INPUT
    %       noiseDataSet    structure which fully describe stochastic process;
    %       count           count of requested samples.
    %
    %   OUTPUT
    %       noise    generated samples.
    %
    narginchk(2, 2);
    
    noise = column_vector_replicate(noiseDataSet.mean, count);
    idxArr = noiseDataSet.idxArr;
    
    for j = 1 : noiseDataSet.N
        ind1 = idxArr(j, 1);
        ind2 = idxArr(j, 2);
        
        switch noiseDataSet.covarianceType
            case 'full'
                sx = chol(noiseDataSet.covariance(ind1:ind2,ind1:ind2))';
            case 'diag'
                sx = diag(sqrt(diag(noiseDataSet.covariance(ind1:ind2, ind1:ind2))));
            case 'sqrt'
                sx = noiseDataSet.covariance(ind1:ind2, ind1:ind2);
            case 'sqrt-diag'
                sx = noiseDataSet.covariance(ind1:ind2, ind1:ind2);
            otherwise
                error('[ sample_combo_gaussian::noiseDataSet ] unknown covarianceType.');
        end
        
        noise(ind1:ind2, :) = noise(ind1:ind2, :) + sx * randn(ind2-ind1+1, count);
    end
end
