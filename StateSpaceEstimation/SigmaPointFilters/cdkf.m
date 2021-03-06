function [ stateNew, stateCovNew, stateNoise, observNoise, internal ] = cdkf(state, covState, stateNoise, observNoise, observation, model, control1, control2 )
    % cdkf. Central Difference Kalman Filter  (Sigma-Point Kalman Filter variant)
    %
    %   [ stateNew, stateCovNew, stateNoise, observNoise, internal ] = cdkf(state, covState, stateNoise, observNoise, observation, model, control1, control2 )
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), u1(k-1)];
    %     z(k) = h[x(k), n(k), u2(k)],
    %
    %   where
    %       x  	is the system state;
    %       v  	the process noise;
    %       n  	the observation noise;
    %       u1 	the exogenous input to the state;
    %       f  	transition function;
    %       u2 	the exogenous input to the state observation function;
    %       z  	the noisy observation of the system.
    %
    %   INPUT
    %         state         state mean at time k-1 ( x(k-1) );
    %         covState     	state covariance at time k-1 ( Px(k-1) );
    %         stateNoise    process noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observNoise   observation noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observation   noisy observations at time k ( z(k) );
    %         model         inference data structure, which fully describes filtration issue (generated by inferenceDataGenerator function);
    %         control1      exogenous input to state transition function at time k-1 ( u1(k-1) );
    %         control2      exogenous input to state observation function at time k ( u2(k) ).
    %
    %   OUTPUT
    %         stateNew          estimates of state starting at time k ( E[x(t)|z(1), z(2), ..., z(t)] for t = k );
    %         stateCovNew       state covariance at time k;
    %         stateNoise        process noise data structure     (possibly updated);
    %         observNoise       observation noise data structure (possibly updated);
    %         internals             <<optional>> internal variables data structure
    %           .meanPredictedState         	predicted state mean ( E[x(t)|z(1), z(2), ..., z(t-1)] for t = k );
    %           .predictedStateCov              predicted state covariance;
    %           .predictedObservMean            predicted observation ( E[z(k)|Z(k-1)] );
    %           .inov                           inovation signal;
    %           .predictedObservCov             inovation covariance;
    %           .filterGain                     filter gain.
    %
    %   Required model fields:
    %         .spkfParams: SPKF parameters = [h] with h: cdkf scale factor / difference step size.
    %
    %% error checking
    if (nargin ~= 8 && nargin ~= 6)
        error('[ cdkf ] Not enough input arguments (should be 6 or 8).');
    end
    
    if (model.stateDimension ~= size(state, 1))
        error('[ cdkf ] Prior state dimension differs from model.stateDimension');
    end
    
    if (model.stateDimension ~= size(covState, 1))
        error('[ cdkf ] Prior state covariance dimension differs from model.stateDimension');
    end
    
    if (model.observationDimension ~= size(observation, 1))
        error('[ cdkf ] Observation dimension differs from model.observationDimension');
    end
    
    %%
    stateDim       = model.stateDimension;
    stateNoiseDim  = model.processNoiseDimension;
    obsNoiseDim    = model.observationNoiseDimension;
    h              = model.spkfParams;
    hh             = h^2;
    
    % sigma-point weights set 1
    w1 = [(hh - stateDim - stateNoiseDim) / hh   1 / (2*hh);
        1 / (2*h)                             sqrt(hh-1) / (2*hh)];
    
    % sigma-point weights set 2
    w2        = w1;
    w2(1, 1)  = (hh - stateDim - obsNoiseDim) / hh;
    
    numSigmaSet1 = 2*(stateDim + stateNoiseDim) + 1;
    numSigmaSet2 = 2*(stateDim + obsNoiseDim) + 1;
    
    if (model.controlInputDimension  == 0)
        control1 = [];
    end
    
    if (model.control2InputDimension == 0)
        control2 = [];
    end
    
    if stateNoise.covariance == 0
        squareRootProcNoiseCov = 0;
    else
        squareRootProcNoiseCov = chol(stateNoise.covariance, 'lower');
    end
    
    squareRootObsNoiseCov  = chol(observNoise.covariance, 'lower');
    squareRootProcCov      = chol(covState, 'lower');
    
    %% time update (prediction)
    if stateNoise.covariance == 0
        sigmaSet1   = cvecrep(state, numSigmaSet1);
        covStateExt = squareRootProcCov;
    else
        sigmaSet1   = cvecrep([state; stateNoise.mean], numSigmaSet1);
        covStateExt = [squareRootProcCov zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) squareRootProcNoiseCov];
    end
    
    sigmaSet1(:, 2 : numSigmaSet1) = sigmaSet1(:, 2:numSigmaSet1) + [h*covStateExt -h*covStateExt];
    
    %% propagate sigma-points through process model
    predictState = model.stateTransitionFun(model, sigmaSet1(1:stateDim, :), sigmaSet1(stateDim + 1 : stateDim + stateNoiseDim, :), control1);
    meanPredictedState = w1(1, 1) * predictState(:, 1) + w1(1, 2) * sum(predictState(:, 2:numSigmaSet1), 2);
    
    a = w1(2, 1) * ( predictState(:, 2:stateDim + stateNoiseDim+1) - predictState(:, stateDim+stateNoiseDim + 2 : numSigmaSet1) ) ;
    b = w1(2, 2) * ( predictState(:, 2:stateDim + stateNoiseDim+1) + predictState(:, stateDim+stateNoiseDim + 2 : numSigmaSet1) - ...
        cvecrep(2*predictState(:, 1), stateDim + stateNoiseDim));
    
    [~, predictStateCov] = qr([a b]', 0);
    predictStateCov= predictStateCov';
    
    sigmaSet2   = cvecrep([meanPredictedState; observNoise.mean], numSigmaSet2);
    crossCovExt = [predictStateCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) squareRootObsNoiseCov];
    sigmaSet2(:, 2:numSigmaSet2) = sigmaSet2(:, 2:numSigmaSet2) + [h*crossCovExt -h*crossCovExt];
    
    %% propagate sigma-points through observation model
    predictObs = model.stateObservationFun(model, sigmaSet2(1:stateDim, :), sigmaSet2(stateDim+1:stateDim+obsNoiseDim, :), control2);
    meanPredictObs = w2(1, 1) * predictObs(:, 1) + w2(1, 2) * sum(predictObs(:, 2:numSigmaSet2), 2);
    c = w2(2,1) * ( predictObs(:, 2:stateDim + obsNoiseDim + 1) - predictObs(:, stateDim+obsNoiseDim + 2:numSigmaSet2) );
    d = w2(2,2) * ( predictObs(:, 2:stateDim + obsNoiseDim + 1) + predictObs(:, stateDim+obsNoiseDim + 2:numSigmaSet2) - cvecrep(2*predictObs(:, 1), stateDim+obsNoiseDim));
    
    % qr([c d]', 0);
    [~, predictObsCov] = qr([c d observNoise.covariance]', 0);
    predictObsCov = predictObsCov';
    
    %% measurement update (correction)
    crossCov = predictStateCov*c(:, 1:stateDim)';
    filterGain = (crossCov / predictObsCov') / predictObsCov;
    
    if isempty(model.innovationModelFunc)
        inov = observation - meanPredictObs;
    else
        inov = model.innovationModelFunc( model, observation, meanPredictObs);
    end
    
    stateNew = meanPredictedState + filterGain*inov;
    
    [~, rCovState] = qr([predictStateCov-filterGain*c(:, 1:stateDim) filterGain*c(:, stateDim+1:end) filterGain*d]', 0);
    rCovState = rCovState';
    stateCovNew = 	rCovState*rCovState';
    
    %% build additional ouptut param (required for debug)
    if nargout > 4
        internal.meanPredictedState    = meanPredictedState;
        internal.predictedStateCov     = predictStateCov;
        internal.predictedObservMean   = meanPredictObs;
        internal.inov                  = inov;
        internal.predictedObservCov    = predictObsCov;
        internal.filterGain            = filterGain;
    end
end