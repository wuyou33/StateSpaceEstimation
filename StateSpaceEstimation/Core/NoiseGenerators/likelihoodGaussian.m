function likelihood = likelihoodGaussian(gaussDataSet, noise, logFlag)
    % likelihoodGaussian. Calculates the likelihood of sample 'noise', given the Gaussian noise model gaussDataSet.
    %
    %   Calculates the likelihood of a 'real world' observation for a given realization or instance of the state variable which has Gaussian noise.
    %   i.e. Calculates the value of P(observation|state).
    %
    %   likelihood = likelihoodGaussian(gaussDataSet, noise, logFlag)
    %
    %   INPUT
    %       gaussDataSet     model of gaussian stochastic process;
    %       noise            values stochastic process;
    %       logFlag     	 <<optional>> determine that log likelihood is used instead of normal likelihood.
    %
    %   OUPUT
    %       likelihood    calculated likelihood.
    %
    narginchk(2, 3);
    
    if nargin == 2
        logFlag = 0;
    end
    
    normFact = (2*pi)^(gaussDataSet.dimension / 2);
    sqrtCov = chol(gaussDataSet.covariance, 'lower');
    xc = sqrtCov \ (noise - cvecrep(gaussDataSet.mean, size(noise, 2)));
    
    if logFlag
        likelihood = -0.5*sum(xc.*xc, 1) - log(normFact*abs(prod(diag(sqrtCov))));
    else
        likelihood = exp(-0.5*sum(xc.*xc, 1))./(normFact*abs(prod(diag(sqrtCov))));
    end
end
