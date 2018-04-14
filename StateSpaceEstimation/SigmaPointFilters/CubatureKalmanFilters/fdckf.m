function [ newState, newCovState, stateNoise, observNoise, internal ] = fdckf( state, covState, stateNoise, observNoise, observation, model, ctrl1, ctrl2)
    % fdckf. Fifth Degree Cubature Kalman Filter. Fifth degree Spherical Simplex-Radial Cubature Kalman Filter (subclass of Sigma Point Kalman Filter).
    % Approximation points (sigma-point) calculated via Spherical Simplex-Radial Rule.
    %
    %  [ newState, newSqrtCov, stateNoise, observNoise, internal ] = fdckf( state, sqrtCovState, stateNoise, observNoise, observation, model, ctrl1, ctrl2)
    %
    %   This filter assumes the following standard state-space model:
    %
    %     x(k) = f[x(k-1), v(k-1), u1(k-1)];
    %     z(k) = h[x(k), n(k), u2(k)],
    %
    %   where
    %       x  - is the system state;
    %       v  - the process noise;
    %       n  - the observation noise;
    %       u1 - the exogenous input to the state;
    %       f  - transition function;
    %       u2 - the exogenous input to the state observation function;
    %       z  - the noisy observation of the system.
    %
    %   INPUT
    %         state         state mean at time k-1 ( x(k-1) );
    %         covState      state covariance at time k-1 ( Px(k-1) );
    %         stateNoise    process noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observNoise   observation noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observation   noisy observations at time k ( z(k) );
    %         model         inference data structure;
    %         ctrl1         exogenous input to state transition function at time k-1 ( u1(k-1) );
    %         ctrl2         exogenous input to state observation function at time k ( u2(k) ).
    %
    %   OUTPUT
    %         newState                  estimates of state starting at time k ( E[x(t)|z(1), z(2), ..., z(t)] for t = k );
    %         newCovState               updated state covariance at time k;
    %         stateNoise                process noise data structure (possibly updated);
    %         observNoise               observation noise data structure (possibly updated);
    %         internal                  <<optional>> internal variables data structure;
    %           .meanPredictedState         	predicted state mean ( E[x(t)|z(1), z(2), ..., z(t-1)] for t = k );
    %           .stateCov                       predicted state covariance matrix at time k;
    %           .predictedObservMean            predicted observation ( E[z(k)|Z(k-1)] );
    %           .inov                           inovation signal;
    %           .observCov                      predicted of Cholesky factor of observation covariance;
    %           .filterGain                     filter gain.
    %
    %% error checking
    if (nargin ~= 8 && nargin ~= 6)
        error('[ fdckf ] Not enough input arguments (should be 6 or 8).');
    end
    
    if (model.stateDimension ~= size(state, 1))
        error('[ fdckf ] Prior state dimension differs from model.stateDimension');
    end
    
    if (model.stateDimension ~= size(covState, 1))
        error('[ fdckf ] Prior state covariance dimension differs from model.stateDimension');
    end
    
    if (model.observationDimension ~= size(observation, 1));
        error('[ fdckf ] Observation dimension differs from model.observationDimension');
    end
    
    %%
    stateDim     = model.stateDimension;
    obsDim       = model.observationDimension;
    numCubPoints = 2*stateDim^2 + 1;
    
    if (model.controlInputDimension == 0)
        ctrl1 = [];
    end
    if (model.control2InputDimension == 0)
        ctrl2 = [];
    end
    
    m_evalFifthDegreeCubatureRule = memoize(@eval_fifth_degree_cubature_rule);
    [points, w] = m_evalFifthDegreeCubatureRule(stateDim);
    w_x = row_vector_replicate(w, stateDim);
    w_z = row_vector_replicate(w, obsDim);
    
    %% calculate cubature points
    offsetPrediction = chol(covState, 'lower');
    cubatureSet  = column_vector_replicate(state, numCubPoints) + offsetPrediction*points;
    
    %% propagate cubature-points through process model
    predictState = model.transition_fun(model, cubatureSet, column_vector_replicate(stateNoise.mean, numCubPoints), ctrl1);
    predictStateMean = predictState*w';
    
    centeredState = (predictState - column_vector_replicate(predictStateMean, numCubPoints));
    predictedStateCov = w_x.*centeredState*centeredState' + stateNoise.covariance;
    
    %% calculate cubature points for measurement
    cubatureSet2 = column_vector_replicate(predictStateMean, numCubPoints) + chol(predictedStateCov, 'lower')*points;    
    
    %% propagate through observation model
    
    predictObs = model.observation_fun(model, cubatureSet2, column_vector_replicate(observNoise.mean, numCubPoints), ctrl2);
    predictObsMean = predictObs*w';
    
    %% measurement update (correction)
    x = (cubatureSet2 - column_vector_replicate(predictStateMean, numCubPoints));
    z = (predictObs - column_vector_replicate(predictObsMean, numCubPoints));
    
    innovationCov = w_z.*z*z'+ observNoise.covariance;
    crossCov = w_x.*x*z';
    filterGain = crossCov*pinv(innovationCov);
    
    if isempty(model.innovation)
        inov = observation - predictObsMean;
    else
        inov = model.innovation( model, observation, predictObsMean);
    end
    
    newState = predictStateMean + filterGain * inov;
    newCovState = predictedStateCov - filterGain*innovationCov*filterGain';
    
    %% build additional ouptut param (required for debug)
    if nargout > 4
        internal.meanPredictedState    = predictStateMean;
        internal.stateCov              = predictedStateCov;
        internal.predictedObservMean   = predictObsMean;
        internal.inov                  = inov;
        internal.observCov             = innovationCov;
        internal.filterGain            = filterGain;
    end
end
