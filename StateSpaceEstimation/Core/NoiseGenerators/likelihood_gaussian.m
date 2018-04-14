function [ likelihood ] = likelihood_gaussian(gaussDataSet, x, logFlag)
    % likelihood_gaussian. Calculates the likelihood of sample 'noise', given the Gaussian noise model gaussDataSet.
    %
    %   Calculates the likelihood of a 'real world' observation for a given realization or instance of the state variable which has Gaussian noise.
    %   i.e. Calculates the value of P(observation|state).
    %
    %   likelihood = likelihood_gaussian(gaussDataSet, noise, logFlag)
    %
    %   INPUT
    %       gaussDataSet     model of gaussian stochastic process;
    %       noise            values stochastic process;
    %       logFlag     	 <<optional>> determine that log likelihood is used instead of normal likelihood.
    %
    %   OUPUT
    %       likelihood    calculated likelihood.
    %%
    narginchk(2, 3);
    
    %%
    if nargin == 2
        logFlag = 0;
    end
    
    switch gaussDataSet.covarianceType
        case {'full', 'diag'}
            sqrtCov = chol(gaussDataSet.covariance, 'lower');
        case {'sqrt', 'sqrt-diag'}
            sqrtCov = gaussDataSet.covariance;
        otherwise
            error('[ likelihood_gaussian::gaussDataSet ] Unknown covariance type ');
    end
    
    arrayLength = size(x, 2);
    
    normfact = (2*pi)^(gaussDataSet.dimension/2);
    xc = x - column_vector_replicate(gaussDataSet.mean, arrayLength);
    foo = sqrtCov \ xc;
    
    if logFlag
        likelihood = -0.5*sum(foo.*foo, 1) - log( normfact*abs(prod(diag(sqrtCov))) );
    else
        likelihood = exp(-0.5*sum(foo.*foo, 1)) ./ ( normfact*abs(prod(diag(sqrtCov))) );
    end
end
