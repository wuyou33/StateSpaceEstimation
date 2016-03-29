function [ llh ] = likelihoodCombo( noiseDataSet, noise, idxVec )
% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.
% 'idxVec' is an optional index vector that can be used to indicate which of the N sub-noise sources should be used.
% to calculate the likelihood... this also requires 'noise' to have the same dimension of the relevant sub-noise source.

    if (nargin == 2)
        idxVec = 1:noiseDataSet.N;
    end

    numNS = length(idxVec);

    [~, nov] = size(noise);
    idxArr = noiseDataSet.idxArr;
    
    llh = ones(1, nov);
    
    for j=1:numNS,

        idx1 = idxArr(idxVec(j),1);
        idx2 = idxArr(idxVec(j),2);
        
        subNoiseDS = noiseDataSet.noiseSources{idxVec(j)};
        llh = llh .* subNoiseDS.likelihood( subNoiseDS, noise(idx1:idx2, :));

    end

    llh = llh + 1e-99; % needed to avoid 0 likelihood (cause ill conditioning)
end