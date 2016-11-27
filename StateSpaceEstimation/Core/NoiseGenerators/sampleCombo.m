function [ noise ] = sampleCombo( noiseDataSet, count )
    % sampleCombo. Generate N (count) samples of a noise source specified by the noiseDataSet data structure (mixture of independent stochastic processes).
    %
    %   [ noise ] = sampleCombo( noiseDataSet, count )
    %
    %   INPUT
    %       noiseDataSet    structure which fully describe stochastic process;
    %       count           count of requested samples.
    %
    %   OUTPUT
    %       noise    generated samples.
    %
    narginchk(2, 2);
    
    noise = zeros(noiseDataSet.dimension, count);
    idxArr = noiseDataSet.idxArr;
    
    for j = 1:noiseDataSet.N
        subNoiseModel = noiseDataSet.noiseSources{j};
        noise(idxArr(j, 1) : idxArr(j, 2), :) = subNoiseModel.sample(subNoiseModel, count);
    end
end
