% Generalized state space model for Doppler navigation system (navigation with using radial velocity relative to the Sun)

%%
function [varargout] = gssmDoppler(func, varargin)
    
    switch func
        
        case 'init'
            varargout{1} = init(varargin{1});
            
        otherwise
            error(['Function ''' func ''' not supported.']);
            
    end
end
%%
function model = init(initArgs)
    model.type                       = 'gssm';                      % object type
    model.tag                        = 'Doppler Navigation System'; % ID tag
    
    model.set_params        = @set_params; % function handle to set_params
    model.transition_fun    = @ffun;        % function handle to state transition function
    model.observation_fun   = @hfun;        % function handle to state observation function
    model.prior             = @prior;       % function handle to the state transition function that calculates P(x(k)|x(k-1)), for particle filter
    model.likelihood        = @likelihood;  % function handle to the observation likelihood function that calculates p(y(k)|x(k)), for particle filter
    model.innovation        = @innovation;  % Function-handle to the innovation model function that calculates the difference between the output
    % of the observation function (hfun) and the actual 'real-world' measurement/observation of that signal,
    % for particle filter
    
    model.stateDimension             = 6;
    model.observationDimension       = 1;
    model.paramDimension             = 0;                              % param
    model.controlInputDimension      = 0;                              % exogenous control input 1 dimension
    model.control2InputDimension     = 0;                              % exogenous control input 2 dimension
    model.processNoiseDimension      = 6;                              % process noise dimension
    model.observationNoiseDimension  = 1;                              % observation noise dimension
    model.params                     = initArgs.initialParams;         % setup parameter vector buffer
    model                            = set_params(model, initArgs.initialParams, initArgs.earthEphemeris, initArgs.sunEphemeris);
    
    % Setup process noise source
    processNoiseArg.type           = 'gaussian';
    processNoiseArg.covarianceType = 'full';
    processNoiseArg.tag            = 'GSSM observation noise observation';
    processNoiseArg.dimension      = model.processNoiseDimension;
    processNoiseArg.mean           = initArgs.stateNoiseMean;
    processNoiseArg.covariance     = initArgs.stateNoiseCovariance;
    
    model.processNoise = generate_noise_model(processNoiseArg);
    
    
    % Setup observation noise source
    observationNoiseArg.type           = 'gaussian';
    observationNoiseArg.covarianceType = 'full';
    observationNoiseArg.tag            = 'GSSM observation noise observation';
    observationNoiseArg.dimension      = model.observationNoiseDimension;
    observationNoiseArg.mean           = initArgs.observationNoiseMean;
    observationNoiseArg.covariance     = initArgs.observationNoiseCovariance;
    
    model.observationNoise = generate_noise_model(observationNoiseArg);
end
%%
function updatedModel = set_params(model, params, earthEphemeris, sunEphemeris)
    % Function to unpack a column vector containing system parameters into specific forms
    % needed by FFUN, HFUN and possibly defined sub-functional objects. Both the vectorized (packed)
    % form of the parameters as well as the unpacked forms are stored within the model data structure.
    
    updatedModel = deep_clone(model);
    updatedModel.params = params;
    
    updatedModel.timeTillCurrentEpoch = params(1);
    updatedModel.sampleTime           = params(2);
    updatedModel.time                 = params(3);
    updatedModel.earthEphemeris       = earthEphemeris;
    updatedModel.sunEphemeris         = sunEphemeris;
end

function newState = ffun(model, state, noise, stateControl)
    % State transition function (system dynamics).
    
    tSpan = [model.time - model.sampleTime; model.time];
    odeFun = @(t,y) equation_of_motion_free_fly(t, y, model.timeTillCurrentEpoch);
    [rn, cn] = size(state);
    
    newState = zeros(rn, cn);
    for i = 1:cn
        [~, tmp]       = ode_euler(odeFun, tSpan, state(:, i), model.sampleTime);
        newState(:, i) = tmp(:, end);
    end
    
    if ~isempty(noise); newState  = newState + noise; end
    
    if ~isempty(stateControl); newState = newState - column_vector_replicate(stateControl, cn); end
end

function observ = hfun(model, state, noise, observationControl)
    % State observation function.
    
    earthEphemeris = model.earthEphemeris;
    sunEphemeris = model.sunEphemeris;
    
    cn = size(state, 2);
    observ = zeros(model.observationDimension, cn);
    
    for i = 1:cn
        observ(:, i) = dopplerShift(state(1:6, i), earthEphemeris, sunEphemeris);
    end
    
    if ~isempty(noise); observ = observ + noise; end
    if ~isempty(observationControl); observ = observ + observationControl; end
end

function tranprior = prior(model, predictedState, state, stateControl, processNoiseDataSet)
    %   Calculates P(nextstate|state).
    
    x = predictedState - ffun(model, state, [], stateControl);
    tranprior = processNoiseDataSet.likelihood( processNoiseDataSet, x);
end

function llh = likelihood(model, observation, state, observationControl, observationNoiseDataSet)
    %   Observation likelihood function.
    
    z = observation - hfun(model, state, [], observationControl);
    llh = observationNoiseDataSet.likelihood( observationNoiseDataSet, z);
end

function innov = innovation(~, observation, predictedObservation) % first argument is a model.
    %   Calculates the innovation signal (difference) between the
    %   output of hfun, i.e. observ (the predicted system observation) and an actual 'real world' observation.
    
    innov = observation - predictedObservation;
end
