function [ newState, newCovState, stateNoise, observNoise, internal ] = cqkf( state, covState, stateNoise, observNoise, observation, model, control1, control2)
    % CQKF Cubature Quadrature Kalman Filter (some kind of High-degree Cubature Kalman Filter)
    %   [ newState, newSvdCovState, stateNoise, observNoise, internal ] = cqkf( state, svdCovState, stateNoise, observNoise, observation, model, control1, control2)
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
    %   Cubature points calculated as intersection of unit hyper-sphere and its axes.
    %   Quadrature points calculated as solution of Chwbychev-Laguerre polynoms with order n' and a = (n / 2 - 1).
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
    if (nargin ~= 8 && nargin ~= 6); error(' [ cqkf ] Not enough input arguments (should be 6 or 8).'); end
    
    if (model.stateDimension ~= size(state, 1)); error('[ cqkf ] Prior state dimension differs from model.stateDimension'); end
    
    if (model.stateDimension ~= size(covState, 1)); error('[ cqkf ] Prior state covariance dimension differs from model.stateDimension'); end
    
    if (model.observationDimension ~= size(observation, 1)); error('[ cqkf ] Observation dimension differs from model.observationDimension'); end
    %%
    stateDim = model.stateDimension;
    obsDim   = model.observationDimension;
    order    = model.cqkfParams(1);
    
    if (model.controlInputDimension == 0); control1 = []; end
    if (model.control2InputDimension == 0); control2 = []; end
    
    [points, w] = cubatureQuadraturePoints(stateDim, order);
               
    %% Evaluate cubature points
    numCubSet1 = 2*stateDim*order;
    weights = elrep(w, stateDim, numCubSet1);
    cubatureSet = cvecrep(state, numCubSet1) + svdDecomposition(covState)*points;
    
    %% Propagate cubature-points through process model
    predictedState = model.stateTransitionFun(model, cubatureSet, cvecrep(stateNoise.mean, numCubSet1), control1);

    predictedStateMean = sum(predictedState.*weights, 2);
    squareRootPredictedStateCov = (predictedState - cvecrep(predictedStateMean, numCubSet1)).*weights;
    predictedStateCov = squareRootPredictedStateCov*squareRootPredictedStateCov' + stateNoise.covariance;
    
    %% Evaluate cubature points for measurement
    weights2 = elrep(w, obsDim, numCubSet1);
    cubatureSet2 = cvecrep(predictedStateMean, numCubSet1) + svdDecomposition(predictedStateCov)*points;
    
    %% Propagate through observation model      
    predictObs = model.stateObservationFun(model, cubatureSet2, cvecrep(observNoise.mean, numCubSet1), control2);
    predictObsMean = sum(predictObs.*weights2, 2);
    
    %% Measurement update
    x = (cubatureSet2 - cvecrep(predictedStateMean, numCubSet1)).*weights;
    z = (predictObs-cvecrep(predictObsMean, numCubSet1)).*weights2;

    innovationCov = z*z'+ observNoise.covariance;
    crossCov = x*z';
    filterGain = crossCov*pinv(innovationCov);

    if isempty(model.innovationModelFunc)
        inov = observation - predictObsMean;
    else
        inov = model.innovationModelFunc(model, observation, predictObsMean);
    end

    newState = predictedStateMean + filterGain*inov;
    newCovState = predictedStateCov - filterGain*innovationCov*filterGain';

    %% additional ouptut param (required for debug)
    internal.meanPredictedState    = predictedStateMean;
    internal.predictedStateCov     = predictedStateCov;
    internal.predictedObservMean   = predictObsMean;
    internal.inov                  = inov;
    internal.predictedObservCov    = innovationCov;
    internal.filterGain            = filterGain;
end
