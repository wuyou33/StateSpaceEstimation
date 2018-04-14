function [ idx ] = resample_mixture_components( xNoiseMixcompsCnt, stateMixcompsCnt, weights, resample_method )
    % resample_mixture_components. Residual resampling for Gaussian Mixture Models.
    %
    %   [ idx ] = resample_mixture_components( mixcompsCnt, weights )
    %
    %   Performs the resampling stage of the SIR GMM algorithm in order (number of samples) steps. Applicable for GSPF algorithm.
    %
    %   INPUT
    %       xNoiseMixcompsCnt   number of state noise gaussian mixture components;
    %       stateMixcompsCnt    number of state gaussian mixture components;
    %       weights             component weigths;
    %       resample_method     reample method.
    %
    %   OUTPUT
    %       idx     resampled indices.
    %
    %%
    mixcompsCnt = xNoiseMixcompsCnt * stateMixcompsCnt;
    resampleIndex = resample(resample_method, weights, mixcompsCnt);
    
    [~, rIdx] = sort(rand(1, mixcompsCnt));    
    rIdx = rIdx(1 : stateMixcompsCnt);
    
    idx  = resampleIndex(rIdx);
end
