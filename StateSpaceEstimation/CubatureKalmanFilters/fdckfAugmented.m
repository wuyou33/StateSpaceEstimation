function [ stateNew, stateCovSvdNew, stateNoise, observNoise, internal ] = fdckfAugmented( state, svdCovState, stateNoise, observNoise, observation, model, ctrl1, ctrl2)
    % fdckfAugmented. Augmented Fifth Degree Cubature Kalman Filter.
    %  [ stateNew, stateCovSvdNew, stateNoise, observNoise, internal ] = fdckfAugmented( state, svdCovState, stateNoise, observNoise, observation, model, ctrl1, ctrl2)
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
    %         state                 state mean at time k-1 ( x(k-1) );
    %         svdCovState           square root factor of matrix through Singular Value Decomposition of state covariance at time k-1;
    %         processNoise          process noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observationNoise      observation noise data structure (must be of type 'gaussian' or 'combo-gaussian');
    %         observation       	noisy observations at time k ( z(k) );
    %         model                 inference data structure, which fully describes filtration issue (generated by inferenceDataGenerator function);
    %         ctrl1                 exogenous input to state transition function at time k-1 ( u1(k-1) );
    %         ctrl2                 exogenous input to state observation function at time k  ( u2(k) ).
    %
    %   OUTPUT
    %         stateNew                  estimates of state at time k ( E[x(t)|z(1), z(2),..., z(t)] for t = k );
    %         stateCovSvdNew            estimate of square root factor of matrix through Singular value decomposition of state covariance at time k;
    %         stateNoise                process noise data structure  (possibly updated);
    %         observNoise               observation noise data structure (possibly updated);
    %         internal                  <<optional>> internal variables data structure;
    %           .meanPredictedState         predicted state mean ( E[x(t)|z(1), z(2), ..., z(t-1)] for t = k );
    %           .sqrtCovState               predicted of Cholesky factor of state covariance matrix at time k;
    %           .predictedObservMean        predicted observation ( E[z(k)|Z(k-1)] )'
    %           .inov                       inovation signal;
    %           .sqrtObservCov              predicted of SVD factor of observation covariance;
    %           .filterGain                 filter gain.
    %
    %% error checking
    if (nargin ~= 8 && nargin ~= 6)
        error('[ fdckfAugmented ] Not enough input arguments (should be 6 or 8).');
    end
    
    if (model.stateDimension ~= size(state, 1))
        error('[ fdckfAugmented ] Prior state dimension differs from model.stateDimension');
    end
    
    if (model.stateDimension ~= size(svdCovState, 1))
        error('[ fdckfAugmented ] Prior state covariance dimension differs from model.stateDimension');
    end
    
    if (model.observationDimension ~= size(observation, 1))
        error('[ fdckfAugmented ] Observation dimension differs from model.observationDimension');
    end
    
    %%
    stateDim      = model.stateDimension;
    obsDim        = model.observationDimension;
    stateNoiseDim = model.processNoiseDimension;
    obsNoiseDim   = model.observationNoiseDimension;
    
    if (model.controlInputDimension == 0)
        ctrl1 = [];
    end
    if (model.control2InputDimension == 0)
        ctrl2 = [];
    end
    
    %% calculate cubature points
    augmentDim   = stateDim + stateNoiseDim;
    numCubPoints = 2*augmentDim^2 + 1;
    
    [cubaturePoints, w] = evalFifthDegreeCubatureRule(augmentDim);
    weights             = rvecrep(w, stateDim);
    
    offset = [svdCovState zeros(stateDim, stateNoiseDim); zeros(stateNoiseDim, stateDim) stateNoise.covariance];
    cubatureSet  = cvecrep([state; stateNoise.mean], numCubPoints) + offset*cubaturePoints;
    
    %% propagate cubature-points through process model
    predictState = model.stateTransitionFun(model, cubatureSet(1:stateDim, :), cubatureSet(stateDim + 1 : stateDim + stateNoiseDim, :), ctrl1);
    predictStateMean = sum(weights.*predictState, 2);
    weightedCenteredSet = (predictState - cvecrep(predictStateMean, numCubPoints)).*weights;
    [~, sqrPredictedCov] = qr(weightedCenteredSet', 0);
    sqrPredictedCov = sqrPredictedCov';
    
    %% calculate cubature points for measurement
    augmentDim = stateDim + obsNoiseDim;
    numCubPointSet2 = 2*augmentDim^2 + 1;
    
    [cubaturePoints, w] = evalFifthDegreeCubatureRule(augmentDim);
    weights2 = rvecrep(w, obsDim);
    offset2 = [sqrPredictedCov zeros(stateDim, obsNoiseDim); zeros(obsNoiseDim, stateDim) observNoise.covariance];
    cubatureSet2 = cvecrep([predictStateMean; observNoise.mean], numCubPointSet2) + offset2*cubaturePoints;
    
    %% propagate through observation model
    predictObs = model.stateObservationFun(model, cubatureSet2(1:stateDim, :), cubatureSet2(stateDim+1:stateDim+obsNoiseDim, :), ctrl2);
    predictObsMean = sum(weights2.*predictObs, 2);
    
    %% measurement update
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
    
    stateNew = predictStateMean + filterGain * inov;
    
    [~, stateCovSvdNew] = qr([(x - filterGain*z)  filterGain*observNoise.covariance]', 0);
    stateCovSvdNew = stateCovSvdNew';
    
    %% build additional ouptut param (required for debug)
    if nargout > 4
        internal.meanPredictedState    = predictState;
        internal.predictedStateCov     = predictStateCov;
        internal.predictedObservMean   = predictedObserv;
        internal.inov                  = inov;
        internal.predictedObservCov    = predictObservCov;
        internal.filterGain            = filterGain;
    end
end
