function llh = likelihoodGamma(noiseDataSet, noise)
%% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.
%%
    llh = noise.^(noiseDataSet.alpha-1) .* exp((-1/noiseDataSet.beta)*noise);
end
