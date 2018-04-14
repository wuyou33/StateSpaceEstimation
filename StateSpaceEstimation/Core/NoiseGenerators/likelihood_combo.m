function [ llh ] = likelihood_combo(noiseDataSet, noise, idxVec)
    % likelihood_combo. Calculates the likelihood of stochastic process (noise).
    %
    % [ llh ] = likelihood_combo(noiseDataSet, noise, idxVec)
    %
    % Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    % 'idxVec' is an optional index vector that can be used to indicate which of the N sub-noise sources should be used to calculate the likelihood.
    % This also requires 'noise' to have the same dimension of the relevant sub-noise source.
    %
    %   INPUT
    %       noiseDataSet    model, which describe stochastic process;
    %       noise           realization (vector of sample) of stochastic process;
    %       idxVec          <<optional>> index vector that can be used to indicate which of the N sub-noise sources should be used to calculate the likelihood.
    %
    %   OUPUT
    %       llh    likelihood vector.
    %
    narginchk(2, 3);
    
    if nargin == 2
        idxVec = 1 : noiseDataSet.N;
    end
    
    noiseSetCount = length(idxVec);
    idxArr = noiseDataSet.idxArr;
    
    llh = ones(1, size(noise, 2));
    
    for j = 1:noiseSetCount
        idx1 = idxArr(idxVec(j), 1);
        idx2 = idxArr(idxVec(j), 2);
        
        subNoiseDS = noiseDataSet.noiseSources{idxVec(j)};
        llh = llh .* subNoiseDS.likelihood(subNoiseDS, noise(idx1:idx2, :));
    end
    
    % needed to avoid 0 likelihood (cause conditioning)
    llh = llh + 1e-99;
end
