function llh = likelihoodGamma(noiseDataSet, noise)
    % likelihoodGamma. Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    %
    %   llh = likelihoodGamma(noiseDataSet, noise)
    %
    %   INPUT
    %       noiseDataSet   inference structure, which fully describe stochastic process;
    %       noise          sample vector of stochastic process.
    %
    %   OUTPUT
    %       llh    likelihood of sample of stochastic process.
    %
    narginchk(2, 2);
    
    llh = noise.^(noiseDataSet.alpha - 1) .* exp((-1/noiseDataSet.beta) * noise);
end
