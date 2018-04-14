function [ noiseDataSet ] = update_combo_gaussian( noiseDataSet )
    % update_combo_gaussian. Updates a 'combination of Gaussian' stochastic processes, which has N Gaussian sub noise sources.
    % The global mean and covariance is updated externally and then this function is called to update the internal sub-noise source structure.
    %
    %   [ noiseDataSet ] = update_combo_gaussian( noiseDataSet )
    %
    %   INPUT
    %       noiseDataSet    structure, which fully describe stochastic process (is a combination of Gaussian noises).
    %
    %   OUTPUT
    %       noiseDataSet    updated structure, which fully describe stochastic process (is a combination of Gaussian noises).
    
    idxArr = noiseDataSet.idxArr;
    
    for j = 1:noiseDataSet.N
        idxRange = idxArr(j, 1) : idxArr(j, 2);
        
        noiseDataSet.noiseSources{j}.mean = noiseDataSet.mean(idxRange, 1);
        noiseDataSet.noiseSources{j}.covariance = noiseDataSet.covariance(idxRange, idxRange);
    end
end
