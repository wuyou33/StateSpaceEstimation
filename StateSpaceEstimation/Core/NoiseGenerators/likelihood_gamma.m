function [ likelihood ] = likelihood_gamma(noiseDataSet, noise)
    % likelihood_gamma. Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    %
    %   likelihood = likelihood_gamma(noiseDataSet, noise)
    %
    %   INPUT
    %       noiseDataSet   inference structure, which fully describe stochastic process;
    %       noise          sample vector of stochastic process.
    %
    %   OUTPUT
    %       likelihood     likelihood of sample of stochastic process.
    %
    narginchk(2, 2);
    
    likelihood = noise.^(noiseDataSet.alpha - 1) .* exp( (-1 / noiseDataSet.beta) * noise );
end
