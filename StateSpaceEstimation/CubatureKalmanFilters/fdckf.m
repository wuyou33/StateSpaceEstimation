function [ newState, newSqrtCov, stateNoise, observNoise, internal ] = fdckf( state, sqrtCovState, stateNoise, observNoise, observation, model, ctrl1, ctrl2)
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
    %         sqrtCovState  square root factor of matrix through Singular Value Decomposition of state covariance at time k-1;
    %         stateNoise    process noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observNoise   observation noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observation   noisy observations at time k ( z(k) );
    %         model         inference data structure;
    %         ctrl1         exogenous input to state transition function at time k-1 ( u1(k-1) );
    %         ctrl2         exogenous input to state observation function at time k ( u2(k) ).
    %
    %   OUTPUT
    %         newState                  estimates of state starting at time k ( E[x(t)|z(1), z(2), ..., z(t)] for t = k );
    %         newSqrtCov                estimate of square root factor of matrix through Singular Value Decomposition of state covariance at time k;
    %         stateNoise                process noise data structure (possibly updated);
    %         observNoise               observation noise data structure (possibly updated);
    %         internal                  <<optional>> internal variables data structure;
    %           .meanPredictedState         	predicted state mean ( E[x(t)|z(1), z(2), ..., z(t-1)] for t = k );
    %           .sqrtCovState                   predicted of Cholesky factor of state covariance matrix at time k;
    %           .predictedObservMean            predicted observation ( E[yz(k)|Z(k-1)] );
    %           .inov                           inovation signal;
    %           .sqrtObservCov                  predicted of Cholesky factor of observation covariance;
    %           .filterGain                     filter gain.
    %
    %% error checking
    if (nargin ~= 8 && nargin ~= 6)
        error('[ fdckf ] Not enough input arguments (should be 6 or 8).');
    end
    
    if (model.stateDimension ~= size(state, 1))
        error('[ fdckf ] Prior state dimension differs from model.stateDimension');
    end
    
    if (model.stateDimension ~= size(sqrtCovState, 1))
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
    
    [cubaturePoints, w] = evalFifthDegreeCubatureRule(stateDim);
    weights = rvecrep(w, stateDim);
    
    %% calculate cubature points
    cubatureSet  = cvecrep(state, numCubPoints) + sqrtCovState*cubaturePoints;
    
    %% propagate cubature-points through process model
    predictState = model.stateTransitionFun(model, cubatureSet, cvecrep(stateNoise.mean, numCubPoints), ctrl1);
    predictStateMean = sum(weights.*predictState, 2);
    weightedCenteredSet = (predictState - cvecrep(predictStateMean, numCubPoints)).*weights;
    [~, sqrPredictedCov] = qr([weightedCenteredSet stateNoise.covariance]', 0);
    sqrPredictedCov = sqrPredictedCov';
    
    %% calculate cubature points for measurement
    cubatureSet2 = cvecrep(predictStateMean, numCubPoints) + sqrPredictedCov*cubaturePoints;
    weights2 = rvecrep(w, obsDim);
    
    %% propagate through observation model
    
    predictObs = model.stateObservationFun(model, cubatureSet2, cvecrep(observNoise.mean, numCubPoints), ctrl2);
    predictObsMean = sum(weights2.*predictObs, 2);
    
    %% measurement update (correction)
    x = (cubatureSet2 - cvecrep(predictStateMean, numCubPoints)) .* weights;
    z = (predictObs - cvecrep(predictObsMean, numCubPoints)) .* weights2;
    
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
    
    [~, qx] = qr([(x - filterGain*z)  filterGain*observNoise.covariance]', 0);
    newSqrtCov = qx';
    
    %% build additional ouptut param (required for debug)
    if nargout > 4
        internal.meanPredictedState    = predictStateMean;
        internal.predictedStateCov     = sqrPredictedCov*sqrPredictedCov';
        internal.predictedObservMean   = predictObsMean;
        internal.inov                  = inov;
        internal.predictedObservCov    = sqrtObsCov*sqrtObsCov';
        internal.filterGain            = filterGain;
    end
end
