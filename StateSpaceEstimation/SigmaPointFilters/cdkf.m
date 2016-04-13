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
    %   Required InferenceDS fields:
    %         .spkfParams: SPKF parameters = [h] with h: CDKF scale factor / difference step size
    
    %% TODO: continue here
    
    
end