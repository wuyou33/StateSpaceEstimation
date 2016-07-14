function [ noise ] = sampleComboGaussian(noiseDataSet, noiseDimension)
%  Generate N samples of a noise source specified by the noiseDataSet data structure

    noise = cvecrep(noiseDataSet.mean, noiseDimension);
    idxArr = noiseDataSet.idxArr;

    for j = 1 : noiseDataSet.N
        ind1 = idxArr(j,1);
        ind2 = idxArr(j,2);
        
        switch noiseDataSet.covarianceType
            case 'full'
                a = chol(noiseDataSet.cov(ind1:ind2,ind1:ind2))';
            case 'diag'
                a = diag(sqrt(diag(noiseDataSet.cov(ind1:ind2, ind1:ind2))));
            case 'sqrt'
                a = noiseDataSet.cov(ind1:ind2, ind1:ind2);
            case 'sqrt-diag'
                a = noiseDataSet.cov(ind1:ind2, ind1:ind2);
            otherwise
                error(' [ sample_gaussian ] unknown cov_type.');
        end
        
        noise(ind1:ind2, :) = noise(ind1:ind2, :) + a * randn(ind2-ind1+1, noiseDimension);
    end
end

