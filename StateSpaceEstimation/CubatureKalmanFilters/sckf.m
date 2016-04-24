function [ newState, newSvdCovState, processNoise, observationNoise, internalVariables ] = sckf( state, svdCovState, processNoise, observationNoise, observation, gssModel, controlProcess, controlObservation )
    % SCKF  Square Root Cubature Kalman Filter (Sigma-Point Kalman Filter variant)
    %
    %   [ newState, newCholCovState, processNoise, observationNoise, internalVariables ] = sckf( state, cholCovState, processNoise, observationNoise, observation, gssModel, controlProcess, controlObservation )
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
    %         svdCovState            square root factor of matrix through Singular value decomposition of state covariance at time k-1
    %         processNoise           process noise data structure     (must be of type 'gaussian' or 'combo-gaussian')
    %         observationNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian')
    %         observation            noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
    %         gssModel               inference data structure
    %         controlProcess         exogenous input to state transition function starting at time k-1 ( u1(k-1),u1(k),...,u1(k+N-2) )
    %         controlObservation     exogenous input to state observation function starting at time k  ( u2(k),u2(k+1),...,u2(k+N-1) )
    %
    %   OUTPUT
    %         newState                   estimates of state starting at time k ( E[x(t)|y(1),y(2),...,y(t)] for t=k,k+1,...,k+N-1 )
    %         newSvdCovState             estimate of square root factor of matrix through Singular value decomposition of state covariance at time k
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
    %% ERROR CHECKING
    if (nargin ~= 8 && nargin ~= 6); error(' [ sckf ] Not enough input arguments (should be 6 or 8).'); end

    if (gssModel.stateDimension ~= size(state, 1)); error('[ sckf ] Prior state dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.stateDimension ~= size(svdCovState, 1)); error('[ sckf ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end

    if (gssModel.observationDimension ~= size(observation, 1)); error('[ sckf ] Observation dimension differs from inferenceDataSet.observationDimension'); end

    %%    
    stateDim         = gssModel.stateDimension;
    procNoiseDim     = gssModel.processNoiseDimension;
    obsNoiseDim      = gssModel.observationNoiseDimension;
    observDim        = gssModel.observationDimension;
        
    augmentDim = stateDim + procNoiseDim;
    numCubPointSet1 = 2*augmentDim;
    
    if (gssModel.controlInputDimension == 0); controlProcess = zeros(0, numCubPointSet1); end
    if (gssModel.control2InputDimension == 0); controlObservation = zeros(0, numCubPointSet1); end
    
    %% Calculate cubature points    
    offsetPrediction = [svdCovState zeros(stateDim, procNoiseDim); zeros(procNoiseDim, stateDim) processNoise.covariance];
    cubatureSet  = cvecrep([state; processNoise.mean], numCubPointSet1) + offsetPrediction*(sqrt(numCubPointSet1/2)*[eye(augmentDim) -eye(augmentDim)]);

    %% Propagate sigma-points through process model
    predictedState = zeros(stateDim, numCubPointSet1);
    for i = 1:numCubPointSet1
        predictedState(:, i) = gssModel.stateTransitionFun(gssModel, cubatureSet(1:stateDim, i), cubatureSet(stateDim+1 : stateDim + procNoiseDim, i), controlProcess(:, i));
    end
    
    predictedStateMean = sum(predictedState, 2) / numCubPointSet1;
    weightedCenteredSet = (predictedState - cvecrep(predictedStateMean, numCubPointSet1)) / sqrt(numCubPointSet1);
    [~, sqrPredictedCov] = qr(weightedCenteredSet', 0);
    sqrPredictedCov = sqrPredictedCov';
    
    %% Calculate cubature points for measurement
    augmentDim = stateDim + obsNoiseDim;
    numCubPointSet2 = 2*augmentDim;
    offsetObs = [sqrPredictedCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) observationNoise.covariance];
    cubatureSet2 = cvecrep([predictedStateMean; observationNoise.mean], numCubPointSet2) + offsetObs*(sqrt(numCubPointSet2/2)*[eye(augmentDim) -eye(augmentDim)]);
    
    %% Propagate through observation model
    predictedObs = zeros(observDim, numCubPointSet2);
    for i = 1:numCubPointSet2
        predictedObs(:, i) = gssModel.stateObservationFun(gssModel, cubatureSet2(1:stateDim, i), cubatureSet2(stateDim+1:stateDim+obsNoiseDim, i), controlObservation(:, i));
    end

    predictedObsMean = sum(predictedObs, 2) / numCubPointSet2;
    
    %% Measurement update
    x = (cubatureSet2(1:stateDim, :)-cvecrep(predictedStateMean, numCubPointSet2)) / sqrt(numCubPointSet2);
    z = (predictedObs-cvecrep(predictedObsMean, numCubPointSet2)) / sqrt(numCubPointSet2);
    
    [~, sqrtObsCov] = qr([z observationNoise.covariance]', 0);
    sqrtObsCov = sqrtObsCov';
    
    crossCov = x*z';
    
    filterGain = (crossCov / sqrtObsCov') / sqrtObsCov;
    
    if isempty(gssModel.innovationModelFunc)
        inov = observation - predictedObsMean;
    else
        inov = gssModel.innovationModelFunc( gssModel, observation, predictedObsMean);
    end
    
    newState = predictedStateMean + filterGain * inov;
    
    [~, newSvdCovState] = qr([(x - filterGain*z)  filterGain*observationNoise.covariance]', 0);
    newSvdCovState = newSvdCovState';
    
    %% additional ouptut param (required for debug)
    internalVariables.meanPredictedState    = predictedStateMean;
    internalVariables.sqrtCovState          = sqrPredictedCov;
    internalVariables.predictedObservMean   = predictedObsMean;
    internalVariables.inov                  = inov;
    internalVariables.sqrtObservCov         = sqrtObsCov;
    internalVariables.filterGain            = filterGain;
end