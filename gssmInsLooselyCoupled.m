% Generalized state space model for loosely coupled inertial navigation
% system & satellite navigation system

%%

function [varargout] = gssmInsLooselyCoupled(func, varargin)
    
    switch func
        
        case 'init'
            model = init(varargin{1});
            varargout{1} = model;
            
        otherwise
            error(['Function ''' func ''' not supported.']);
            
    end
end
%%
function model = init(initArgs)
    
    model.type                       = 'gssm';                         % object type
    model.tag                        = 'Loosely Coupled NS';           % ID tag
    
    model.setParams                  = @setparams;                     % function handle to SETPARAMS
    model.stateTransitionFun         = @ffun;                          % function handle to state transition function
    model.stateObservationFun        = @hfun;                          % function handle to state observation function
    
    model.stateTransitionPriorFun    = @prior;                         % function handle to the state transition function that calculates P(x(k)|x(k-1)),
    % is a requirement for particle filter
    
    model.observationLikelihoodFun   = @likelihood;                    % function handle to the observation likelihood function that calculates p(y(k)|x(k)),
    % is a requirement for particle filter
    
    model.innovationModelFunc        = @innovation;                    % Function-handle to the innovation model function that calculates the difference between the output
    % of the observation function (hfun) and the actual 'real-world' measurement/observation of that signal,
    % is a requirement for particle filter
    
    model.stateDimension             = 22;
    model.observationDimension       = 6;
    model.paramDimension             = 22;                             % param count = 7, but each param is 3D vector or quaternion
    model.controlInputDimension      = 10;                             % exogenous control input 1 dimension
    model.control2InputDimension     = 0;                              % exogenous control input 2 dimension
    model.processNoiseDimension      = 22;                             % process noise dimension
    model.observationNoiseDimension  = 6;                              % observation noise dimension
    model.params                     = initArgs.initialParams;         %  setup parameter vector buffer
    model                            = setparams(model, initArgs.initialParams);
    
    % Setup process noise source
    
    processNoiseArg.type           = 'gaussian';
    processNoiseArg.covarianceType = 'full';
    processNoiseArg.tag            = 'GSSM process noise source';
    processNoiseArg.dimension      = model.processNoiseDimension;
    processNoiseArg.mean           = initArgs.processNoiseMean;
    processNoiseArg.covariance     = initArgs.processNoiseCovariance;
    
    model.processNoise = generateNoiseDataSet(processNoiseArg);
    
    % Setup observation noise source
    
    observationNoiseArg.type           = 'gaussian';
    observationNoiseArg.covarianceType = 'full';
    observationNoiseArg.tag            = 'GSSM observation noise observation';
    observationNoiseArg.dimension      = model.observationNoiseDimension;
    observationNoiseArg.mean           = initArgs.observationNoiseMean;
    observationNoiseArg.covariance     = initArgs.observationNoiseCovariance;
    
    model.observationNoise = generateNoiseDataSet(observationNoiseArg);
end
%%
function model = setparams(model, params, idxVector)
    
    % Function to unpack a column vector containing system parameters into specific forms
    % needed by FFUN, HFUN and possibly defined sub-functional objects. Both the vectorized (packed)
    % form of the parameters as well as the unpacked forms are stored within the model data structure.
    % idxVector is an optional argument which indicates which parameters should be updated. This can
    % be used to only modify a subset of the total system parameters.
    
    switch nargin
        
        case 2 % update all
            model.params = params;
        case 3 % update subset of params
            model.params(idxVector) = params;
        otherwise
            error('[ setparams ] Incorrect number of input arguments.');
    end
    
    model.accelerationBiasMu(1:3, 1)        = model.params(1:3);
    model.accelerationBiasSigma(1:3, 1)     = model.params(4:6);
    model.gyroBiasMu(1:3, 1)                = model.params(7:9);
    model.gyroBiasSigma(1:3, 1)             = model.params(10:12);
    model.acceleration(1:3, 1)              = model.params(13:15);
    model.angularVelocity(1:3, 1)           = model.params(16:18);
    model.quaternion(1:4, 1)                = model.params(19:22);
    model.sampleTime                        = model.params(23);
    model.time                              = model.params(24);
end

function newState = ffun(model, state, noise, stateControl)
    % State transition function (system dynamics).
    a = model.acceleration;
    w = model.angularVelocity;
    q = model.quaternion;
    tSpan = [model.time - model.sampleTime; model.time];
    odeFun = @(t, y) SinsDynamicEquation(t, y, a, w, q);
    [rn, cn] = size(state);
    
    newState = zeros(rn, cn);
    for i = 1:cn
        %     [~, tmp] = ode45(odeFun, tSpan, state, odeset('MaxStep', model.sampleTime));
        %     newState = tmp(end, :)';    
        [~, tmp]       = odeEuler(odeFun, tSpan, state(:, i), model.sampleTime);
        newState(:, i) = tmp(:, end);
    end
    
    if ~isempty(noise)
        newState  = newState + cvecrep([ones(10, 1); sqrt(model.sampleTime)*ones(6, 1); ones(6, 1)], cn).*noise;
%         newState(11:16) = newState(11:16) + sqrt(model.sampleTime)*([model.accelerationBiasSigma; model.gyroBiasSigma].*randn(6, 1));
    end
    
    newState(7:10, :) = quaternionNormalize(newState(7:10, :));
    
    if ~isempty(stateControl)
        newState(1:6, :)  = newState(1:6, :) - cvecrep(stateControl(1:6), cn);
        newState(7:10, :) = quaternionNormalize(  quaternionMultiply( newState(7:10, :), quaternionConj(cvecrep(stateControl(7:10), cn), 2) )  );
    end
end

function observ = hfun(~, state, noise, observationControl) % first argument is a model.
    % State observation function.
    observ = state(1:6, :);
    
    if ~isempty(noise)
        observ = observ + noise;
    end
    
    if ~isempty(observationControl)
        observ = observ + cvecrep(observationControl, size(state, 2));
    end
end

function pr = prior(model, predictedState, state, stateControl, processNoiseDataSet)
    % Calculates P(nextstate|state).    
    x = predictedState - ffun(model, state, [], stateControl);        
    pr = processNoiseDataSet.likelihood(processNoiseDataSet, x);
end

function llh = likelihood(model, observation, state, observationControl, observationNoiseDataSet)
    %   Observation likelihood function
    z = observation - hfun(model, state, [], observationControl);    
    llh = observationNoiseDataSet.likelihood(observationNoiseDataSet, z);
end

function innov = innovation(~, observation, predictedObservation) % first argument is a model.
    %   Calculates the innovation signal (difference) between the
    %   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual 'real world' observation.
    innov = observation - predictedObservation;
end