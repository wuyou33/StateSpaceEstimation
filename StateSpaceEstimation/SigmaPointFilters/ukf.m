function [newState, newCovState, processNoise, observationNoise, internalVariables] = ukf(state, covState, processNoise, observationNoise, observation, inferenceDataSet, controlProcess, controlObservation)
% UKF  Unscented Kalman Filter
%
%   [newState, newCovState, processNoise, observationNoise, internalVariables] = ukf(state, covState, processNoise, observationNoise, observation, controlProcess, controlObservation, inferenceDataSet)
%
%   This filter assumes the following standard state-space model:
%
%     x(k) = f[x(k-1), v(k-1), U1(k-1)]
%     y(k) = h[x(k), n(k), U2(k)]
%
%   where: 
%       x is the system state, 
%       v the process noise, 
%       n the observation noise, 
%       u1 the exogenous input to the state
%       f the transition function, 
%       u2 the exogenous input to the state observation function 
%       y the noisy observation of the system.
%
%   INPUT
%         state                  state mean at time k-1          ( xh(k-1) )
%         covState               state covariance at time k-1    ( Px(k-1) )
%         processNoise           process noise data structure     (must be of type 'gaussian' or 'combo-gaussian')
%         observationNoise       observation noise data structure (must be of type 'gaussian' or 'combo-gaussian')
%         observation            noisy observations starting at time k ( y(k),y(k+1),...,y(k+N-1) )
%         controlProcess         exogenous input to state transition function starting at time k-1 ( u1(k-1),u1(k),...,u1(k+N-2) )
%         controlObservation     exogenous input to state observation function starting at time k  ( u2(k),u2(k+1),...,u2(k+N-1) )
%         inferenceDataSet       inference data structure generated by GENINFDS function.
%
%   OUTPUT
%         newState               estimates of state starting at time k ( E[x(t)|y(1),y(2),...,y(t)] for t=k,k+1,...,k+N-1 )
%         newState               state covariance
%         processNoise           process noise data structure     (possibly updated)
%         observationNoise       observation noise data structure (possibly updated)
%
%         internalVariables      <<optional>> internal variables data structure
%           .xh_                    predicted state mean ( E[x(t)|y(1),y(2),..y(t-1)] for t=k,k+1,...,k+N-1 )
%           .Px_                    predicted state covariance
%           .yh_                    predicted observation ( E[y(k)|Y(k-1)] )
%           .inov                   inovation signal
%           .Pinov                  inovation covariance
%           .KG                     Kalman gain
%
%   Required InferenceDS fields:
%         .spkfParams            SPKF parameters = [alpha beta kappa] with
%                                   alpha  :  UKF scale factor
%                                   beta   :  UKF covariance correction factor
%                                   kappa  :  UKF secondary scaling parameter
%%

stateDim         = inferenceDataSet.stateDimension;                   % extract state dimension
procNoiseDim     = inferenceDataSet.processNoiseDimension;            % extract process noise dimension
obserNoiseDim    = inferenceDataSet.observationNoiseDimension;        % extract observation noise dimension

%% ERROR CHECKING
if (nargin ~= 8 && nargin ~= 6); error(' [ ukf ] Not enough input arguments (should be 6 or 8).'); end

if (inferenceDataSet.stateDimension ~= size(state, 1)); error('[ ukf ] Prior state dimension differs from inferenceDataSet.stateDimension'); end

if (inferenceDataSet.stateDimension ~= size(covState, 1)); error('[ ukf ] Prior state covariance dimension differs from inferenceDataSet.stateDimension'); end

if (inferenceDataSet.observationDimension ~= size(observation, 1)); error('[ ukf ] Observation dimension differs from inferenceDataSet.observationDimension'); end

%% Get UKF scaling parameters
alpha = inferenceDataSet.spkfParams(1);
beta  = inferenceDataSet.spkfParams(2);
kappa = inferenceDataSet.spkfParams(3);

augmentStateDim = stateDim + procNoiseDim + obserNoiseDim;
numSigmaPoints  = 2*augmentStateDim + 1;
kappa           = alpha^2*(augmentStateDim + kappa)-augmentStateDim;

if (inferenceDataSet.controlInputDimension == 0); 
    ctrlProc = zeros(0, numSigmaPoints);
else
    ctrlProc = cvecrep(controlProcess, numSigmaPoints);
end

if (inferenceDataSet.control2InputDimension == 0); 
    ctrlObserv = zeros(0, numSigmaPoints);
else
    ctrlObserv = cvecrep(controlObservation, numSigmaPoints);
end

weights    = [kappa 0.5 0] / (augmentStateDim + kappa);
weights(3) = weights(1) + (1 - alpha^2) + beta;

%% generate sigma point set
sigmaPoints  = cvecrep([state; processNoise.mean; observationNoise.mean], numSigmaPoints);
covStateExt  = [chol(covState)' zeros(stateDim, procNoiseDim); zeros(procNoiseDim, stateDim) chol(processNoise.covariance)'];
offset       = sqrt(augmentStateDim + kappa) * ([covStateExt zeros(stateDim+procNoiseDim, obserNoiseDim); zeros(obserNoiseDim, stateDim + procNoiseDim) chol(observationNoise.covariance)']);
fullOffset   = [offset -offset];
sigmaPoints(:, 2 : numSigmaPoints) = sigmaPoints(:, 2:numSigmaPoints) + fullOffset;

%% propagate sigma-points through process model
predictedState = inferenceDataSet.stateTransitionFun(inferenceDataSet, sigmaPoints(1:stateDim, :), sigmaPoints(stateDim+1 : stateDim + procNoiseDim, :), ctrlProc);
meanPredictedState = weights(1) * predictedState(:, 1) + weights(2)*sum(predictedState(:, 2:numSigmaPoints), 2);
centredPredictedState = predictedState - cvecrep(meanPredictedState, numSigmaPoints);
predictedStateCov = weights(3)*centredPredictedState(:, 1)*centredPredictedState(:, 1)' + weights(2)*centredPredictedState(:, 2:numSigmaPoints)*centredPredictedState(:, 2:numSigmaPoints)';

%% propagate through observation model
predictedObserv = inferenceDataSet.stateObservationFun(inferenceDataSet, predictedState, sigmaPoints(stateDim+procNoiseDim+1 : stateDim+procNoiseDim+obserNoiseDim, :), ctrlObserv);
predictedObservMean = weights(1)*predictedObserv(:, 1) + weights(2)*sum(predictedObserv(:, 2:numSigmaPoints), 2);
centredPredicatedObserv = predictedObserv - cvecrep(predictedObservMean, numSigmaPoints);
predictedObservCov  = weights(3)*centredPredicatedObserv(:, 1)*centredPredicatedObserv(:, 1)' + weights(2)*centredPredicatedObserv(:, 2:numSigmaPoints)*centredPredicatedObserv(:, 2:numSigmaPoints)';

%% measurement update
crossCov = weights(3)*centredPredictedState(:,1)*centredPredicatedObserv(:,1)' + weights(2)*centredPredictedState(:,2:numSigmaPoints)*centredPredicatedObserv(:,2:numSigmaPoints)';
filterGain = crossCov / predictedObservCov;

if isempty(inferenceDataSet.innovation)
    inov = observation - predictedObservMean;
else
    inov = inferenceDataSet.innovationState(inferenceDataSet, observation, predictedObservMean);
end

newState = meanPredictedState + filterGain*inov;
newCovState = predictedStateCov - filterGain*predictedObservCov*filterGain';

%% additional ouptut param (required for debug)
internalVariables.meanPredictedState    = meanPredictedState;
internalVariables.predictedStateCov     = predictedStateCov;
internalVariables.predictedObservMean   = predictedObservMean;
internalVariables.inov                  = inov;
internalVariables.predictedObservCov    = predictedObservCov;
internalVariables.filterGain            = filterGain;