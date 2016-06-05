% Generalized state space model for loosely coupled inertial navigation
% system & satellite navigation system

%%

function [varargout] = gssmInsSns(func, varargin)

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
    model.tag                        = 'Loosely Coupled INS & SNS';    % ID tag

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
    model.controlInputDimension      = 0;                              % exogenous control input 1 dimension
    model.control2InputDimension     = 0;                              % exogenous control input 2 dimension
    model.processNoiseDimension      = 6;                              % process noise dimension
    model.observationNoiseDimension  = 6;                              % observation noise dimension
    model.params                     = initArgs.initialParams;         %  setup parameter vector buffer
    model                            = setparams(model, ...
                                            initArgs.initialParams);

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

    model.accelerationBiasMu     = model.params(1:3);
    model.accelerationBiasSigma  = model.params(4:6);
    model.gyroBiasMu             = model.params(7:9);
    model.gyroBiasSigma          = model.params(10:12);
    model.acceleration           = model.params(13:15);
    model.angularVelocity        = model.params(16:18);
    model.quaternion             = model.params(19:22);
    model.sampleTime             = model.params(23);
    model.time                   = model.params(24);
end

function newState = ffun(model, state, noise, stateControl)
    % State transition function (system dynamics).    
    newState = rungeKuttaFourOrderWithFixedStep(...
        @(t, y) SinsDynamicEquation(t, ...
            y, ...
            model.acceleration', ...
            model.angularVelocity', ...
            model.quaternion ), ...
        state', ...
        model.time, ...
        model.sampleTime);

    newState(7:10) = quaternionNormalize(newState(7:10));
    
    if ~isempty(noise)
        newState(17:22) = newState(17:22) + sqrt(model.sampleTime)*noise';
    end
    
    if ~isempty(stateControl)
        newState = newState + stateControl;
    end
end

function observ = hfun(~, state, noise, observationControl) % first argument is a model.
    % State observation function.
    observ = state(1:6);

    if ~isempty(noise)
        observ = observ + noise(1,:);
    end
    
    if ~isempty(observationControl)
        observ = observ + observationControl;
    end
end

function tranprior = prior(model, predictedState, state, stateControl, processNoiseDataSet)
%   Calculates P(nextstate|state).
    x = predictedState - ffun(model, state, [], stateControl);    

    tranprior = processNoiseDataSet.likelihood( processNoiseDataSet, x);
end

function llh = likelihood(model, observation, state, observationControl, observationNoiseDataSet)
%   Observation likelihood function
    x = observation - hfun(model, state, [], observationControl);

    llh = observationNoiseDataSet.likelihood( observationNoiseDataSet, x);
end

function innov = innovation(~, observation, predictedObservation) % first argument is a model.
%   Calculates the innovation signal (difference) between the
%   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual 'real world' observation OBS.
    innov = observation - predictedObservation;
end