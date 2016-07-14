function [ newState, newSvdCovState, stateNoise, observNoise, internal ] = cqkf( state, svdCovState, stateNoise, observNoise, observation, model, control1, control2)
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
    
    if (model.stateDimension ~= size(state, 1)); error('[ cqkf ] Prior state dimension differs from inferenceDataSet.stateDimension'); end
    
    if (model.stateDimension ~= size(svdCovState, 1)); error('[ cqkf ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end
    
    if (model.observationDimension ~= size(observation, 1)); error('[ cqkf ] Observation dimension differs from inferenceDataSet.observationDimension'); end
    %%
    stateDim = model.stateDimension;
    obsDim   = model.observationDimension;
    order    = model.cqkfParams(1);
        
    if (model.controlInputDimension == 0); control1 = []; end
    if (model.control2InputDimension == 0); control2 = []; end
    
    cubaturePoints   = intersectUnitVectorHyperSphere(stateDim);
    [ quadraturePoints, weigths ] = laguerreQuadratureRule(order, alpha);
    points = zeros(2*stateDim*order, stateDim);
    for i = 1:order
        points(stateDim*(1-1) + 1, stateDim*i, :) = sqrt(2)*cubaturePoints*quadraturePoints(i);
    end
    
    
end
