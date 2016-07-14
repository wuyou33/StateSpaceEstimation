function [ newState, newCovState, processNoise, observationNoise, internalVariables ] = cdkf( state, covState, processNoise, observationNoise, observation, gssModel, controlProcess, controlObservation )
    % [xh, Px, pNoise, oNoise, InternalVariablesDS] = cdkf(state, Pstate, pNoise, oNoise, obs, U1, U2, InferenceDS)
    % CDKF:  Central Difference Kalman Filter  (Sigma-Point Kalman Filter variant)
    %
    %   [ newState, newCovState, processNoise, observationNoise, internalVariables ] = cdkf( state, covState, processNoise, observationNoise, observation, gssModel, controlProcess, controlObservation )
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), U1(k-1)]
    %     y(k) = h[x(k), n(k), U2(k)]
    %
    %   where
    %       x  - is the system state,
    %       v  - the process noise,
    %       n  - the observation noise,
    %       U1 - the exogenous input to the state
    %       f  - transition function,
    %       U2 - the exogenous input to the state observation function
    %       y  - the noisy observation of the system.
    %
    %   INPUT
    %         state                  state mean at time k-1          ( xh(k-1) )
    %         covState               state covariance at time k-1    ( Px(k-1) )
    %         processNoise           process noise data structure     (must be of type 'gaussian' or 'combo-gaussian')
    %         observationNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian')
    %         observation            noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
    %         gssModel               inference data structure
    %         controlProcess         exogenous input to state transition function starting at time k-1 ( u1(k-1),u1(k),...,u1(k+N-2) )
    %         controlObservation     exogenous input to state observation function starting at time k  ( u2(k),u2(k+1),...,u2(k+N-1) )
    %
    %   OUTPUT
    %         newState                   estimates of state starting at time k ( E[x(t)|y(1),y(2),...,y(t)] for t=k,k+1,...,k+N-1 )
    %         newCovState                state covariance
    %         processNoise               process noise data structure     (possibly updated)
    %         observationNoise           observation noise data structure (possibly updated)
    %         internalVariables             <<optional>> internal variables data structure
    %           .meanPredictedState         	predicted state mean ( E[x(t)|y(1),y(2),..y(t-1)] for t=k,k+1,...,k+N-1 )
    %           .predictedStateCov              predicted state covariance
    %           .predictedObservMean            predicted observation ( E[y(k)|Y(k-1)] )
    %           .inov                           inovation signal
    %           .predictedObservCov             inovation covariance
    %           .filterGain                     Kalman gain
    %
    %   Required gssModel fields:
    %         .spkfParams: SPKF parameters = [h] with h: CDKF scale factor / difference step size
    
    %% ERROR CHECKING
    if (nargin ~= 8 && nargin ~= 6); error(' [ cdkf ] Not enough input arguments (should be 6 or 8).'); end

    if (gssModel.stateDimension ~= size(state, 1)); error('[ cdkf ] Prior state dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.stateDimension ~= size(covState, 1)); error('[ cdkf ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.observationDimension ~= size(observation, 1)); error('[ cdkf ] Observation dimension differs from inferenceDataSet.observationDimension'); end
    
    %%        
    stateDim       = gssModel.stateDimension;    
    stateNoiseDim  = gssModel.processNoiseDimension;
    obsNoiseDim    = gssModel.observationNoiseDimension;    
    h              = gssModel.spkfParams;
    hh             = h^2;
    
    % sigma-point weights set 1
    w1 = [(hh - stateDim - stateNoiseDim) / hh   1 / (2*hh);
           1 / (2*h)                             sqrt(hh-1) / (2*hh)];
    
    % sigma-point weights set 2
    w2       = w1;
    w2(1, 1) = (hh - stateDim - obsNoiseDim) / hh;
        
    numSigmaSet1 = 2*(stateDim + stateNoiseDim) + 1;
    numSigmaSet2 = 2*(stateDim + obsNoiseDim) + 1;
    
    if (gssModel.controlInputDimension  == 0); controlProcess = []; end
    
    if (gssModel.control2InputDimension == 0); controlObservation = []; end
    
    if processNoise.covariance == 0
        squareRootProcNoiseCov = 0;
    else
        squareRootProcNoiseCov = chol(processNoise.covariance, 'lower');
    end
    
    squareRootObsNoiseCov  = chol(observationNoise.covariance, 'lower');
    squareRootProcCov      = chol(covState, 'lower');
    
    %% time update
    if processNoise.covariance == 0
        sigmaSet1 = cvecrep(state, numSigmaSet1);
        covStateExt    = squareRootProcCov;
    else
        sigmaSet1   = cvecrep([state; processNoise.mean], numSigmaSet1);
        covStateExt = [squareRootProcCov zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) squareRootProcNoiseCov];
    end    
    
    sigmaSet1(:, 2 : numSigmaSet1) = sigmaSet1(:, 2:numSigmaSet1) + [h*covStateExt -h*covStateExt];
    
    %% propagate sigma-points through process model
    predictState = gssModel.stateTransitionFun(gssModel, sigmaSet1(1:stateDim, :), sigmaSet1(stateDim + 1 : stateDim + stateNoiseDim, :), controlProcess);               
    meanPredictedState = w1(1, 1) * predictState(:, 1) + w1(1, 2) * sum(predictState(:, 2:numSigmaSet1), 2);
    
    a = w1(2, 1) * ( predictState(:, 2:stateDim + stateNoiseDim+1) - predictState(:, stateDim+stateNoiseDim + 2 : numSigmaSet1) ) ;
    b = w1(2, 2) * ( predictState(:, 2:stateDim + stateNoiseDim+1) + predictState(:, stateDim+stateNoiseDim + 2 : numSigmaSet1) - ...
        cvecrep(2*predictState(:, 1), stateDim + stateNoiseDim));
    
    [~, predictStateCov] = qr([a b]', 0);
    predictStateCov= predictStateCov';

    sigmaSet2   = cvecrep([meanPredictedState; observationNoise.mean], numSigmaSet2);   
    crossCovExt = [predictStateCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) squareRootObsNoiseCov];    
    sigmaSet2(:, 2:numSigmaSet2) = sigmaSet2(:, 2:numSigmaSet2) + [h*crossCovExt -h*crossCovExt];
    
    %% propagate sigma-points through observation model
    predictObs = gssModel.stateObservationFun(gssModel, sigmaSet2(1:stateDim, :), sigmaSet2(stateDim+1:stateDim+obsNoiseDim, :), controlObservation);    
    meanPredictObs = w2(1, 1) * predictObs(:, 1) + w2(1, 2) * sum(predictObs(:, 2:numSigmaSet2), 2);
    c = w2(2,1) * ( predictObs(:, 2:stateDim + obsNoiseDim + 1) - predictObs(:, stateDim+obsNoiseDim + 2:numSigmaSet2) );
    d = w2(2,2) * ( predictObs(:, 2:stateDim + obsNoiseDim + 1) + predictObs(:, stateDim+obsNoiseDim + 2:numSigmaSet2) - cvecrep(2*predictObs(:, 1), stateDim+obsNoiseDim));
    
    [~, predictObsCov] = qr([c d observationNoise.covariance]', 0); % qr([c d]', 0);
    predictObsCov = predictObsCov';
    
    %% measurement update    
    crossCov = predictStateCov*c(:, 1:stateDim)';
    filterGain = (crossCov / predictObsCov') / predictObsCov;

    if isempty(gssModel.innovationModelFunc)
        inov = observation - meanPredictObs;
    else
        inov = gssModel.innovationModelFunc( gssModel, observation, meanPredictObs);
    end

    newState = meanPredictedState + filterGain*inov;

    [~, rCovState] = qr([predictStateCov-filterGain*c(:, 1:stateDim) filterGain*c(:, stateDim+1:end) filterGain*d]', 0);
    rCovState = rCovState';
    newCovState = 	rCovState*rCovState';
    
    %% additional ouptut param (required for debug)
    internalVariables.meanPredictedState    = meanPredictedState;
    internalVariables.predictedStateCov     = predictStateCov;
    internalVariables.predictedObservMean   = meanPredictObs;
    internalVariables.inov                  = inov;
    internalVariables.predictedObservCov    = predictObsCov;
    internalVariables.filterGain            = filterGain;
end