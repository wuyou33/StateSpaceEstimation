function [ noise ] = updateComboGaussian( noise )
% Updates a 'combination Gaussian' noise source which has N Gaussian sub noise sources. The global mean and covariance
% is updated externally and then this function is called to update the internal sub-noise source structure.%

    idxArr = noise.idxArr;

    for j=1:noise.N,
        ind1 = idxArr(j,1);
        ind2 = idxArr(j,2);
        idxRange = ind1:ind2;
        
        noise.noiseSources{j}.mean = noise.mean(idxRange, 1);
        noise.noiseSources{j}.covariance = noise.covariance(idxRange, idxRange);
    end
end