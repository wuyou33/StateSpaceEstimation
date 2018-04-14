function [ likelihood ] = likelihood_gmm( noiseDataSet, noise )
    % likelihood_gmm. Calculates the likelihood of sample 'noise', given the noise model noiseDataSet.
    %
    %   Calculates the likelihood of a 'real world' observation for a given realization or instance of the state variable,
    %   which describes as mixture of Gaussian stochastic processes, i.e. Calculates the value of P(observation|state).
    %
    %   [ likelihood ] = likelihood_gmm( noiseDataSet, noise )
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
    [ ~, ~, likelihood] = gmm_probability(noiseDataSet, noise);
end
