function likelihood = likelihoodGaussian(gausDataSet, x, logFlag)
    % likelihoodGaussian
    %   likelihood = likelihoodGaussian(gausDataSet, x, logFlag)
    %   INPUT:
    %       gausDataSet - model of gaussian noise;
    %       x           - values of state vector;
    %       logFlag     - determine that log likelihood is used instead of normal likelihood.
    %   OUPUT:
    %       likelihood  - calculated likelihood.
    %   Calculates the likelihood of a 'real world' observation 'OBS' for a given
    %   realization or instance of the state variable STATE
    %   which has Gaussian noise. i.e. Calculates the value of P(OBS|STATE).
    %
    %%
        
    if nargin ~=2 && nargin ~= 3
        error(' [ likelihoodGaussian ] Not enough inputs.');
    elseif nargin == 2
        logFlag = 0;
    end
        
    normFact = (2*pi)^(gausDataSet.dimension / 2);
    sqrtCov = chol(gausDataSet.covariance, 'lower');
    xc = sqrtCov \ (x - cvecrep(gausDataSet.mean, size(x, 2)));
    
    if logFlag
        likelihood = -0.5*sum(xc.*xc, 1) - log(normFact*abs(prod(diag(sqrtCov))));
    else
        likelihood = exp(-0.5*sum(xc.*xc, 1))./(normFact*abs(prod(diag(sqrtCov))));
    end
end
