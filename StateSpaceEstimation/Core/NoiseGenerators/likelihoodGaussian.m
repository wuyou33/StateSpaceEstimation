function likelihood = likelihoodGaussian(gausDataSet, x, logFlag)
%{
Calculates the likelihood of a 'real world' observation 'OBS' for a given 
realization or instance of the state variable STATE 
which has Gaussian noise. i.e. Calculates the value of P(OBS|STATE). 
%}
%%
    switch nargin
        case 3
        case 2
            logFlag = 0;
        otherwise
            error(' [ likelihoodGaussian ] Not enough inputs.');
    end

    normFact = (2*pi)^(gausDataSet.dimension/2);
    covDecomposition = chol(gausDataSet.covariance)';
    foo = covDecomposition \ (x - cvecrep(gausDataSet.mean, size(x, 2)));
    if logFlag
        likelihood = -0.5*sum(foo.*foo, 1) - log(normFact*abs(prod(diag(covDecomposition))));
    else
        likelihood = exp(-0.5*sum(foo.*foo, 1))./(normFact*abs(prod(diag(covDecomposition))));
    end
end
