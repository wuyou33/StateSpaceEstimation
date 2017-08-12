function [ idx ] = resampleMixtureComponents( xNoiseMixcompsCnt, stateMixcompsCnt, weights )
    % resampleMixtureComponents. Residual resampling for Gaussian Mixture Models.
    %
    %   [ idx ] = resampleMixtureComponents( mixcompsCnt, weights )
    %
    %   Performs the resampling stage of the SIR GMM algorithm in order (number of samples) steps. Applicable for GSPF algorithm.
    %
    %   INPUT
    %       xNoiseMixcompsCnt   number of state noise gaussian mixture components;
    %       stateMixcompsCnt    number of state gaussian mixture components;
    %       weights             component weigths.
    %
    %   OUTPUT
    %       idx     resampled indices.
    %
    %%
    mixcompsCnt = xNoiseMixcompsCnt * stateMixcompsCnt;
    resampleIndex = residualResample(1:mixcompsCnt, weights);
    [~, rIdx] = sort(rand(1, mixcompsCnt));
    rIdx = rIdx(1 : stateMixcompsCnt);
    idx  = resampleIndex(rIdx);
end
