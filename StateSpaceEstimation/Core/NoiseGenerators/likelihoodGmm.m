function [ llh ] = likelihoodGmm( noiseDataSet, noise )
% Calculates the likelihood of sample 'noise', given the noise model NoiseDS.

    [ ~, ~, llh] = gmmProbability(noiseDataSet, noise); % prior, likelihood, evidence
end
