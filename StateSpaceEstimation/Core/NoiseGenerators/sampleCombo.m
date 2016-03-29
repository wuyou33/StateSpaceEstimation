function [ noise ] = sampleCombo( noiseDataSet, noiseDimension )
%  Generate N samples of a noise source specified by the NoiseDS data structure
    noise = zeros(noiseDataSet.dimension, noiseDimension);
    idxArr = noiseDataSet.idxArr;

    for j = 1:noiseDataSet.N
        subNoiseDS = noiseDataSet.noiseSources{j};
        noise(idxArr(j,1) : idxArr(j,2), :) = subNoiseDS.sample(subNoiseDS, noiseDimension);
    end
end