function [estimate, dataSet, stateNoise, observNoise] = pf(dataSet, stateNoise, observNoise, observation, control1, control2, model)
    % PF  Generic Particle Filter
    %   This filter is also known as the 'Bootstrap Particle Filter' or the 'Condensation Algorithm'
    %
    %   [estimate, dataSet, stateNoise, observNoise] = pf(dataSet, stateNoise, observNoise, observation, control1, control2, model)
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
    %
    %%
    if nargin ~= 7; error(' [ pf ] Incorrect number of input arguments.'); end
    
    stateDim    = model.stateDimension;
    num         = dataSet.particlesNum;
    particles   = dataSet.particles;
    weights     = dataSet.weights;
    threshold   = round(model.resampleThreshold * num);
    normWeights = cvecrep(1/num, num);
    
    if (model.controlInputDimension == 0); control1 = []; end
    
    if (model.control2InputDimension == 0); control2 = []; end
    
    xNoise = stateNoise.sample(stateNoise, num);
    %% propagate particles
    particlesPred = model.stateTransitionFun(model, particles, xNoise, control1);
    
    %% evaluate importance weights
    likelihood = model.likelihoodStateFun(model, cvecrep(observation, num), particlesPred, control2, observNoise) + 1e-99;
    weights = weights .* likelihood;
    weights = weights / sum(weights);
    
    %% resample
    effectiveSetSize = 1 / sum(weights.^2); % calculate effective particle set size
    
    if (effectiveSetSize < threshold) % resample if effectiveSetSize is below threshold
        outIndex  = residualResample(1:num, weights);
        particles = particlesPred(:, outIndex);
        weights   = normWeights;
    else
        particles  = particlesPred;
    end
    
    %% caculate estimate
    if strcmp(model.estimateType, 'mean')
        estimate = sum(rvecrep(weights, stateDim).*particles, 2);
    elseif strcmp(model.estimateType, 'median')
        estimate = median(rvecrep(weights, stateDim).*particles, 2);
    else
        error(' [ pf ] Unknown estimate type.');
    end
    
    dataSet.particles = particles;
    dataSet.weights   = weights;
end
