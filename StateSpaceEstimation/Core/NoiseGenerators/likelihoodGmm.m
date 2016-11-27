function [ likelihood ] = likelihoodGmm( noiseDataSet, noise )
    % likelihoodGmm. Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    %
    %   Calculates the likelihood of a 'real world' observation for a given realization or instance of the state variable,
    %   which describes as mixture of Gaussian stochastic processes, i.e. Calculates the value of P(observation|state).
    %
    %   [ likelihood ] = likelihoodGmm( noiseDataSet, noise )
    %
    %   INPUT
    %       noiseDataSet     model of gaussian stochastic process;
    %       noise            values stochastic process;
    %
    %   OUPUT
    %       likelihood    calculated likelihood.
    %
    narginchk(2, 2);
    
    % prior, likelihood, evidence
    [ ~, ~, likelihood] = gmmProbability(noiseDataSet, noise);
end
