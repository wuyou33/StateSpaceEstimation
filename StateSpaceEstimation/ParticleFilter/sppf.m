function [estimate, dataSet, stateNoise, observNoise] = sppf(dataSet, stateNoise, observNoise, observation, control1, control2, model)
    % sppf  Sigma-Point Particle Filter.
    %   This hybrid particle filter uses a sigma-point Kalman filter (SRUKF, SRCDKF) or Cubature Kalman filter (SCKF) for proposal distribution generation
    %   and is an extension of the original "Unscented Particle Filter.
    %
    %   [estimate, dataSet, processNoise, observationNoise] = sppf(dataSet, processNoise, observationNoise, observation, controlProc, controlObs, model)
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), U1(k-1)]
    %     z(k) = h[x(k), n(k), U2(k)]
    %
    %   where:
    %       x is the system state,
    %       v the process noise,
    %       n the observation noise,
    %       u1 the exogenous input to the state
    %       f the transition function,
    %       u2 the exogenous input to the state observation function
    %       y the noisy observation of the system.
    %
    %   INPUT
    %         dataSet           particle filter data structure. Contains set of particles as well as their corresponding weights.
    %         stateNoise        process noise data structure
    %         observNoise       observation noise data structure
    %         observation       noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
    %         control1          exogenous input to state transition function starting at time k-1 ( u1(k-1),u1(k),...,u1(k+N-2) )
    %         control2          exogenous input to state observation function starting at time k  ( u2(k),u2(k+1),...,u2(k+N-1) )
    %         gssModel          inference data structure.
    %
    %   OUTPUT
    %         estimate          state estimate generated from posterior distribution of state given all observation. Type of
    %                           estimate is specified by 'InferenceDS.estimateType'
    %         dataSet           updated Particle filter data structure. Contains set of particles as well as their corresponding weights.
    %         stateNoise        process noise data structure     (possibly updated)
    %         observNoise       observation noise data structure (possibly updated)
    %
    %   dataSet fields:
    %         .particlesNum        (scalar) number of particles
    %         .particles           (statedim-by-N matrix) particle buffer
    %         .weights             (1-by-N r-vector) particle weights
    %
    %   Required gssModel fields:
    %         .estimateType        (string) Estimate type : 'mean', 'mode', etc.
    %         .resampleThreshold   (scalar) If the ratio of the 'effective particle set size' to the total number of particles
    %                                       drop below this threshold  i.e.  (nEfective / particlesNum) < resampleThreshold
    %                                       the particles will be resampled.  (nEfective is always less than or equal to particlesNum)
    %%
    if nargin ~= 7; error('[ sppf ] Not enough input arguments.'); end
    
    %%
    stateDim        = model.stateDimension;
    num             = dataSet.particlesNum;
    particles       = dataSet.particles;
    sqrtCov         = dataSet.particlesCov;
    procNoiseKF     = dataSet.processNoise;
    obsNoiseKF      = dataSet.observationNoise;
    weights         = dataSet.weights;
    threshold       = round(num*model.resampleThreshold);
    
    normWeights = cvecrep(1 / num, num);
    
    if (model.controlInputDimension == 0); control1 = []; end    
    if (model.control2InputDimension == 0); control2 = []; end
    
    sqrtCovPred = zeros(stateDim, stateDim, num);
    stateNew    = zeros(stateDim, num);
    statePred   = zeros(stateDim, num);
    proposal    = zeros(1, num);
    normFact    = (2*pi)^(stateDim / 2);
    
    %% time update
    
    switch model.spkfType
        case 'srukf'
            predict = @(x, s, xNoise, zNoise) srukf(x, s, xNoise, zNoise, observation, model, control1, control2);
        case 'sckf'
            predict = @(x, s, xNoise, zNoise) sckf(x, s, xNoise, zNoise, observation, model, control1, control2);
        case 'srcdkf'
            predict = @(x, s, xNoise, zNoise) srcdkf(x, s, xNoise, zNoise, observation, model, control1, control2);
        otherwise
            error('[ sppf ] Unknown inner filter type.');
    end
    
    for i = 1:num
        [stateNew(:, i), sqrtCovPred(:, :, i), procNoiseKF, obsNoiseKF] = predict(particles(:, i), sqrtCov(:, :, i), procNoiseKF, obsNoiseKF);
        statePred(:, i) = stateNew(:, i) + sqrtCovPred(:, :, i)*randn(stateDim, 1);
    end
    
    %% evalutate importance weights
    
    prior = model.stateTransitionPriorFun(model, statePred, particles, control1, stateNoise) + 1e-99;
    
    likelihood = model.likelihoodStateFun(model, cvecrep(observation, num), statePred, control2, observNoise) + 1e-99;
    
    difState = statePred - stateNew;
    
    for i=1:num
        covFact = sqrtCovPred(:, :, i);
        expData = covFact \ difState(:, i);
        proposal(i) = exp(-0.5*(expData'*expData)) / abs( normFact * prod( diag(covFact) ) ) + 1e-99;
        weights(i) = weights(i) * likelihood(i) * prior(i) / proposal(i);
    end
    
    weights = weights / sum(weights);
    
    if strcmp(model.estimateType, 'mean')
        estimate = sum( weights(ones(1, stateDim), :) .* statePred, 2);
    else
        error('[ sppf ] Unknown estimate type.');
    end
    
    %% resample
    
    effectiveSize = 1 / sum(weights.^2);
    
    if (effectiveSize < threshold)
        outIndex  = residualResample(1:num,weights);
        particles = statePred(:, outIndex);
        
        for i = 1:num
            sqrtCov(:, :, i) = sqrtCovPred(:, :, outIndex(i));
        end
        
        weights = normWeights;
    else
        particles   = statePred;
        sqrtCov     = sqrtCovPred;
    end
    
    dataSet.particles           = particles;
    dataSet.particlesCov        = sqrtCov;
    dataSet.weights             = weights;
    dataSet.processNoise        = procNoiseKF;
    dataSet.observationNoise    = obsNoiseKF;
    
end
