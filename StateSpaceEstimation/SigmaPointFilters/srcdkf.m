function [ newState, newCholCovState, processNoise, observationNoise, internalVariables ] = srcdkf( state, cholCovState, processNoise, observationNoise, observation, gssModel, controlProcess, controlObservation )
    % SRCDKF  Square Root Central Difference Kalman Filter (Sigma-Point Kalman Filter variant)
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
    %         cholCovState           upper triangle metrix from Cholesky decomposition of state covariance at time k-1    ( chol(Px(k-1)) )
    %         processNoise           process noise data structure     (must be of type 'gaussian' or 'combo-gaussian')
    %         observationNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian')
    %         observation            noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
    %         gssModel               inference data structure
    %         controlProcess         exogenous input to state transition function starting at time k-1 ( u1(k-1),u1(k),...,u1(k+N-2) )
    %         controlObservation     exogenous input to state observation function starting at time k  ( u2(k),u2(k+1),...,u2(k+N-1) )
    %
    %   OUTPUT
    %         newState                   estimates of state starting at time k ( E[x(t)|y(1),y(2),...,y(t)] for t=k,k+1,...,k+N-1 )
    %         newCholCovState            estimate of Cholesky factor of state covariance matrix at time k
    %         processNoise               process noise data structure     (possibly updated)
    %         observationNoise           observation noise data structure (possibly updated)
    %         internalVariables             <<optional>> internal variables data structure
    %           .meanPredictedState         	predicted state mean ( E[x(t)|y(1),y(2),..y(t-1)] for t=k,k+1,...,k+N-1 )
    %           .sqrtCovState                   predicted of Cholesky factor of state covariance matrix at time k
    %           .predictedObservMean            predicted observation ( E[y(k)|Y(k-1)] )
    %           .inov                           inovation signal
    %           .sqrtObservCov                  predicted of Cholesky factor of observation covariance
    %           .filterGain                     Kalman gain
    %
    %   Required gssModel fields:
    %         .spkfParams: SPKF parameters = [h] with h: CDKF scale factor / difference step size
    
    %% ERROR CHECKING
    if (nargin ~= 8 && nargin ~= 6); error(' [ srcdkf ] Not enough input arguments (should be 6 or 8).'); end

    if (gssModel.stateDimension ~= size(state, 1)); error('[ srcdkf ] Prior state dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.stateDimension ~= size(cholCovState, 1)); error('[ srcdkf ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.observationDimension ~= size(observation, 1)); error('[ srcdkf ] Observation dimension differs from inferenceDataSet.observationDimension'); end
    
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
    w2      = w1;
    w2(1,1) = (hh - stateDim - obsNoiseDim)/hh;
        
    numSigmaPointSet1 = 2*(stateDim + stateNoiseDim) + 1;
    numSigmaPointSet2 = 2*(stateDim + obsNoiseDim) + 1;
    
    if (gssModel.controlInputDimension  == 0); controlProcess = []; end
    
    if (gssModel.control2InputDimension == 0); controlObservation = zeros(0, numSigmaPointSet2); end
    
    %% time update
    if processNoise.covariance == 0
        sigmaPointSet1 = cvecrep(state, numSigmaPointSet1);    
        offsetSet1     = cholCovState;
    else
        sigmaPointSet1 = cvecrep([state; processNoise.mean], numSigmaPointSet1);    
        offsetSet1     = [cholCovState zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) processNoise.covariance];    
    end
    sigmaPointSet1(:, 2:numSigmaPointSet1) = sigmaPointSet1(:, 2:numSigmaPointSet1) + [h*offsetSet1 -h*offsetSet1];    
    
    %% propagate sigma-points through process model
    predictState = gssModel.stateTransitionFun(gssModel, sigmaPointSet1(1:stateDim, :), sigmaPointSet1(stateDim + 1 : stateDim + stateNoiseDim, :), controlProcess);    
    meanPredictState = w1(1, 1) * predictState(:, 1) + w1(1, 2) * sum(predictState(:, 2:numSigmaPointSet1), 2);
    
    a = w1(2, 1) * ( predictState(:, 2:stateDim + stateNoiseDim+1) - predictState(:, stateDim+stateNoiseDim + 2 : numSigmaPointSet1) ) ;
    b = w1(2, 2) * ( predictState(:, 2:stateDim + stateNoiseDim+1) + predictState(:, stateDim+stateNoiseDim + 2 : numSigmaPointSet1) - ...
        cvecrep(2*predictState(:, 1), stateDim + stateNoiseDim));
    
    [~, sqrtPredictStateCov] = qr([a b]', 0);
    sqrtPredictStateCov= sqrtPredictStateCov';
    
    sigmaPointSet2 = cvecrep([meanPredictState; observationNoise.mean], numSigmaPointSet2);       
    offsetSet2     = [sqrtPredictStateCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) observationNoise.covariance];    
    sigmaPointSet2(:, 2:numSigmaPointSet2) = sigmaPointSet2(:,2:numSigmaPointSet2) + [h*offsetSet2 -h*offsetSet2];
    
    %% propagate sigma-points through observation model
    predictedObs = gssModel.stateObservationFun(gssModel, sigmaPointSet2(1:stateDim, :), sigmaPointSet2(stateDim+1:stateDim+obsNoiseDim, :), controlObservation);        
    meanPredictedObs = w2(1, 1) * predictedObs(:, 1) + w2(1, 2) * sum(predictedObs(:, 2:numSigmaPointSet2), 2);
    c = w2(2,1) * ( predictedObs(:, 2:stateDim + obsNoiseDim + 1) - predictedObs(:, stateDim + obsNoiseDim + 2 : numSigmaPointSet2) );
    d = w2(2,2) * ( predictedObs(:, 2:stateDim + obsNoiseDim + 1) + predictedObs(:, stateDim + obsNoiseDim + 2 : numSigmaPointSet2) - cvecrep(2*predictedObs(:, 1), stateDim + obsNoiseDim));
    
    [~, predictedObservCov] = qr([c d observationNoise.covariance]', 0); % qr([c d]', 0);
    predictedObservCov = predictedObservCov';
    
    %% measurement update    
    crossCov = sqrtPredictStateCov*c(:, 1:stateDim)';
    filterGain = (crossCov / predictedObservCov') / predictedObservCov;

    if isempty(gssModel.innovationModelFunc)
        inov = observation - meanPredictedObs;
    else
        inov = gssModel.innovationModelFunc( gssModel, observation, meanPredictedObs);
    end

    newState = meanPredictState + filterGain*inov;

    [~, newCholCovState] = qr([sqrtPredictStateCov-filterGain*c(:, 1:stateDim) filterGain*c(:, stateDim+1:end) filterGain*d]', 0);
    newCholCovState = newCholCovState';
    
    %% additional ouptut param (required for debug)
    internalVariables.meanPredictedState    = meanPredictState;
    internalVariables.predictedStateCov     = sqrtPredictStateCov;
    internalVariables.predictedObservMean   = meanPredictedObs;
    internalVariables.inov                  = inov;
    internalVariables.predictedObservCov    = predictedObservCov;
    internalVariables.filterGain            = filterGain;
end