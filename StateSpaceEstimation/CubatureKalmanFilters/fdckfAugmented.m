function [ newState, newSvdCovState, stateNoise, observNoise, internal ] = fdckfAugmented( state, svdCovState, stateNoise, observNoise, observation, model, control1, control2)
    % FDCKF Augmented Fifth Degree Cubature Kalman Filter
    %  [ newState, newSvdCovState, stateNoise, observNoise, internal ] = sckf( state, svdCovState, stateNoise, observNoise, observation, gssModel, control1, control2)
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), u1(k-1)]
    %     y(k) = h[x(k), n(k), u2(k)]
    %
    %   where
    %       x  - is the system state,
    %       v  - the process noise,
    %       n  - the observation noise,
    %       u1 - the exogenous input to the state
    %       f  - transition function,
    %       u2 - the exogenous input to the state observation function
    %       y  - the noisy observation of the system.
    %
    %   INPUT
    %         state                  state mean at time k-1          ( xh(k-1) )
    %         svdCovState            square root factor of matrix through Singular value decomposition of state covariance at time k-1
    %         processNoise           process noise data structure     (must be of type 'gaussian' or 'combo-gaussian')
    %         observationNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian')
    %         observation            noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
    %         model                  inference data structure
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
    if (nargin ~= 8 && nargin ~= 6); error(' [ fdckfAugmented ] Not enough input arguments (should be 6 or 8).'); end
    
    if (model.stateDimension ~= size(state, 1)); error('[ fdckfAugmented ] Prior state dimension differs from inferenceDataSet.stateDimension'); end
    
    if (model.stateDimension ~= size(svdCovState, 1)); error('[ fdckfAugmented ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end
    
    if (model.observationDimension ~= size(observation, 1)); error('[ fdckfAugmented ] Observation dimension differs from inferenceDataSet.observationDimension'); end
    
    %%
    stateDim      = model.stateDimension;
    obsDim        = model.observationDimension;
    stateNoiseDim = model.processNoiseDimension;
    obsNoiseDim   = model.observationNoiseDimension;
    
    if (model.controlInputDimension == 0); control1 = []; end
    if (model.control2InputDimension == 0); control2 = []; end
    
    %% Calculate cubature points
    augmentDim   = stateDim + stateNoiseDim;
    numCubPoints = 2*augmentDim^2 + 1;
    
    [cubaturePoints, w] = evalFifthDegreeCubatureRule(augmentDim);
    weights      = rvecrep(w, stateDim);
    
    offset = [svdCovState zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) stateNoise.covariance];
    cubatureSet  = cvecrep([state; stateNoise.mean], numCubPoints) + offset*cubaturePoints;
    
    %% Propagate cubature-points through process model
    predictState = model.stateTransitionFun(model, cubatureSet(1:stateDim, :), cubatureSet(stateDim + 1 : stateDim + stateNoiseDim, :), control1);
    predictStateMean = sum(weights.*predictState, 2);
    weightedCenteredSet = (predictState - cvecrep(predictStateMean, numCubPoints)).*weights;
    [~, sqrPredictedCov] = qr(weightedCenteredSet', 0);
    sqrPredictedCov = sqrPredictedCov';
    
    %% Calculate cubature points for measurement
    augmentDim = stateDim + obsNoiseDim;
    numCubPointSet2 = 2*augmentDim^2 + 1;
    
    [cubaturePoints, w] = evalFifthDegreeCubatureRule(augmentDim);
    weights2 = rvecrep(w, obsDim);
    offset2 = [sqrPredictedCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) observNoise.covariance];
    cubatureSet2 = cvecrep([predictStateMean; observNoise.mean], numCubPointSet2) + offset2*cubaturePoints;
    
    %% Propagate through observation model
    predictObs = model.stateObservationFun(model, cubatureSet2(1:stateDim, :), cubatureSet2(stateDim+1:stateDim+obsNoiseDim, :), control2);
    predictObsMean = sum(weights2.*predictObs, 2);
    
    %% Measurement update
    x = (cubatureSet2(1:stateDim, :)-cvecrep(predictStateMean, numCubPointSet2)).*rvecrep(w, stateDim);
    z = (predictObs-cvecrep(predictObsMean, numCubPointSet2)).*weights2;
    
    [~, sqrtObsCov] = qr([z observNoise.covariance]', 0);
    sqrtObsCov = sqrtObsCov';
    crossCov = x*z';
    
    filterGain = (crossCov / sqrtObsCov') / sqrtObsCov;
    
    if isempty(model.innovationModelFunc)
        inov = observation - predictObsMean;
    else
        inov = model.innovationModelFunc( model, observation, predictObsMean);
    end
    
    newState = predictStateMean + filterGain * inov;
    
    [~, newSvdCovState] = qr([(x - filterGain*z)  filterGain*observNoise.covariance]', 0);
    newSvdCovState = newSvdCovState';
    
    %% additional ouptut param (required for debug)
    internal.meanPredictedState    = predictStateMean;
    internal.sqrtCovState          = sqrPredictedCov;
    internal.predictedObservMean   = predictObsMean;
    internal.inov                  = inov;
    internal.sqrtObservCov         = sqrtObsCov;
    internal.filterGain            = filterGain;
end
