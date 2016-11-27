function [stateNew, stateCovNew, stateNoise, observNoise, internals] = kf(state, covState, stateNoise, observNoise, observation, model, control1, control2)
    % kf. Kalman Filter (Linear Kalman Filter).
    %
    %   [stateNew, stateCovNew, stateNoise, observNoise, internals] = kf(state, covState, stateNoise, observNoise, observation, model, control1, control2)
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), u1(k-1)];
    %     z(k) = h[x(k), n(k), u2(k)],
    %
    %   where:
    %       x  - is the system state;
    %       v  - the process noise;
    %       n  - the observation noise;
    %       u1 - the exogenous input to the state;
    %       f  - the transition function;
    %       u2 - the exogenous input to the state observation function;
    %       z  - the noisy observation of the system.
    %
    %   INPUT:
    %         state             state mean at time k-1 ( x(k-1) );
    %         covState          state covariance at time k-1 ( Px(k-1) );
    %         stateNoise        process noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observation       noisy observations starting at time k ( z(k) );
    %         model             inference data structure, which fully describes filtration issue;
    %         control1          exogenous input to state transition function starting at time k-1 ( u1(k-1) );
    %         control2          exogenous input to state observation function starting at time k ( u2(k) ).
    %
    %   OUTPUT
    %         stateNew          estimates of state starting at time k ( E[x(t)|z(1), z(2), ..., z(t)] for t = k );
    %         stateCovNew       state covariance;
    %         stateNoise        process noise data structure (possibly updated);
    %         observNoise       observation noise data structure (possibly updated);
    %         internals             <<optional>> internal variables data structure
    %           .meanPredictedState         	predicted state mean ( E[x(t)|z(1), z(2), ..., z(t-1)] for t = k );
    %           .predictedStateCov              predicted state covariance;
    %           .predictedObservMean            predicted observation ( E[z(k)|Z(k-1)] );
    %           .inov                           inovation signal;
    %           .predictedObservCov             inovation covariance;
    %           .filterGain                     filter gain.
    %
    %% error checking
    
    if (nargin ~= 8 && nargin ~= 6)
        error('[ kf ]  Not enough input arguments (should be 6 or 8).');
    end
    
    if (model.stateDimension ~= size(state, 1))
        error('[ kf ] Prior state dimension differs from model.stateDimension');
    end
    
    if ( model.stateDimension ~= size(covState, 1))
        error('[ kf ] Prior state covariance dimension differs from model.stateDimension');
    end
    
    if (model.observationDimension ~= size(observation, 1))
        error('[ kf ] Observation dimension differs from model.observationDimension');
    end
    %%
    if (model.controlInputDimension == 0)
        control1 = [];
    end
    
    if (model.control2InputDimension == 0);
        control2 = [];
    end
    
    %% time update (prediction)
    F = model.linearize(model, state, stateNoise.mean, [], control1, [], 'F');
    G = model.linearize(model, state, stateNoise.mean, [], control1, [], 'G');
    
    statePredict    = model.stateTransitionFun(model, state, stateNoise.mean, control1);
    covStatePredict = F * covState * F' + G * stateNoise.covariance * G';
    
    %% measurement update (correction)
    C = model.linearize(model, statePredict, [], observNoise.mean, [], control2, 'C');
    H = model.linearize(model, statePredict, [], observNoise.mean, [], control2, 'H');
    
    covObsPredict = C * covStatePredict * C' + H * observNoise.covariance * H';
    filterGain    = covStatePredict * C' / covObsPredict;
    
    observPredict  = model.stateObservationFun(model, statePredict, observNoise.mean, control2);
    
    if isempty(model.innovation)
        inov = observation - observPredict;
    else
        inov = model.innovationModelFunc(model, observation, observPredict);
    end
    
    stateNew    = statePredict + filterGain * inov;
    stateCovNew = covStatePredict - filterGain * covObsPredict * filterGain';
    
    %% build additional ouptut param (required for debug)
    if nargout > 4
        internals.meanPredictedState    = statePredict;
        internals.predictedStateCov     = covStatePredict;
        internals.predictedObservMean   = observPredict;
        internals.inov                  = inov;
        internals.predictedObservCov    = covObsPredict;
        internals.filterGain            = filterGain;
    end
end
