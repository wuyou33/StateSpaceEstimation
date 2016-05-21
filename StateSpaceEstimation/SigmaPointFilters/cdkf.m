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
    weights1 = [(hh - stateDim - stateNoiseDim)/hh   1/(2*hh);
                1/(2*h)                     sqrt(hh-1)/(2*hh)];
    
    % sigma-point weights set 2
    weights2      = weights1;
    weights2(1,1) = (hh - stateDim - obsNoiseDim)/hh;
        
    numSigmaPointSet1 = 2*(stateDim + stateNoiseDim) + 1;
    numSigmaPointSet2 = 2*(stateDim + obsNoiseDim) + 1;
    
    if (gssModel.controlInputDimension  == 0); controlProcess = zeros(0, numSigmaPointSet1); end
    if (gssModel.control2InputDimension == 0); controlObservation = zeros(0, numSigmaPointSet2); end
    
    if processNoise.covariance == 0
        squareRootProcNoiseCov = 0;
    else
        squareRootProcNoiseCov = chol(processNoise.covariance)';
    end
    squareRootObsNoiseCov  = chol(observationNoise.covariance)';
    squareRootProcCov      = chol(covState)';
    
    %% time update
    if processNoise.covariance == 0
        sigmaPointSet1 = cvecrep(state, numSigmaPointSet1);
        covStateExt    = squareRootProcCov;
    else
        sigmaPointSet1 = cvecrep([state; processNoise.mean], numSigmaPointSet1);
        covStateExt    = [squareRootProcCov zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) squareRootProcNoiseCov];
    end    
    
    sigmaPointSet1(:,2:numSigmaPointSet1) = sigmaPointSet1(:,2:numSigmaPointSet1) + [h*covStateExt -h*covStateExt];
    
    %% propagate sigma-points through process model
    predictedState = zeros(stateDim, numSigmaPointSet1);
    for i = 1:numSigmaPointSet1                               
        predictedState(:, i) = gssModel.stateTransitionFun(gssModel, sigmaPointSet1(1:stateDim, i), sigmaPointSet1(stateDim + 1 : stateDim + stateNoiseDim, i), controlProcess(:, i));
    end
    
    meanPredictedState = weights1(1, 1) * predictedState(:, 1) + weights1(1, 2) * sum(predictedState(:, 2:numSigmaPointSet1), 2);
    
    a = weights1(2, 1) * ( predictedState(:, 2:stateDim + stateNoiseDim+1) - predictedState(:, stateDim+stateNoiseDim + 2 : numSigmaPointSet1) ) ;
    b = weights1(2, 2) * ( predictedState(:, 2:stateDim + stateNoiseDim+1) + predictedState(:, stateDim+stateNoiseDim + 2 : numSigmaPointSet1) - ...
        cvecrep(2*predictedState(:, 1), stateDim + stateNoiseDim));
    
    [~, predictedStateCov] = qr([a b]', 0);
    predictedStateCov= predictedStateCov';

    sigmaPointSet2 = cvecrep([meanPredictedState; observationNoise.mean], numSigmaPointSet2);   
    crossCovExt    = [predictedStateCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) squareRootObsNoiseCov];    
    sigmaPointSet2(:,2:numSigmaPointSet2) = sigmaPointSet2(:,2:numSigmaPointSet2) + [h*crossCovExt -h*crossCovExt];
    
    %% propagate sigma-points through observation model
    predictedObs = zeros(obsNoiseDim, numSigmaPointSet2);
    for i = 1:numSigmaPointSet2                                
        predictedObs(:, i) = gssModel.stateObservationFun(gssModel, sigmaPointSet2(1:stateDim, i), sigmaPointSet2(stateDim+1:stateDim+obsNoiseDim, i), controlObservation(:, i));
    end
    
    meanPredictedObs = weights2(1, 1) * predictedObs(:, 1) + weights2(1, 2) * sum(predictedObs(:, 2:numSigmaPointSet2), 2);
    c = weights2(2,1) * ( predictedObs(:, 2:stateDim+obsNoiseDim + 1) - predictedObs(:, stateDim+obsNoiseDim + 2:numSigmaPointSet2) );
    d = weights2(2,2) * ( predictedObs(:, 2:stateDim+obsNoiseDim + 1) + predictedObs(:, stateDim+obsNoiseDim + 2:numSigmaPointSet2) - cvecrep(2*predictedObs(:, 1), stateDim+obsNoiseDim));
    
    [~, predictedObservCov] = qr([c d observationNoise.covariance]', 0); % qr([c d]', 0);
    predictedObservCov = predictedObservCov';
    
    %% measurement update    
    crossCov = predictedStateCov*c(:, 1:stateDim)';
    filterGain = (crossCov / predictedObservCov') / predictedObservCov;

    if isempty(gssModel.innovationModelFunc)
        inov = observation - meanPredictedObs;
    else
        inov = gssModel.innovationModelFunc( gssModel, observation, meanPredictedObs);
    end

    newState = meanPredictedState + filterGain*inov;

    [~, rCovState] = qr([predictedStateCov-filterGain*c(:, 1:stateDim) filterGain*c(:, stateDim+1:end) filterGain*d]', 0);
    rCovState = rCovState';
    newCovState = 	rCovState*rCovState';
    
    %% additional ouptut param (required for debug)
    internalVariables.meanPredictedState    = meanPredictedState;
    internalVariables.predictedStateCov     = predictedStateCov;
    internalVariables.predictedObservMean   = meanPredictedObs;
    internalVariables.inov                  = inov;
    internalVariables.predictedObservCov    = predictedObservCov;
    internalVariables.filterGain            = filterGain;
end