% Generalized state space model for loosely coupled inertial navigation
% system & satellite navigation system

%%

function [varargout] = gssmXNav(func, varargin)

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
    model.tag                        = 'X-Ray Navigation System';      % ID tag

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

    model.stateDimension             = 6;
    model.observationDimension       = 7;
    model.paramDimension             = 0;                             % param count = 7, but each param is 3D vector or quaternion
    model.controlInputDimension      = 0;                              % exogenous control input 1 dimension
    model.control2InputDimension     = 0;                              % exogenous control input 2 dimension
    model.processNoiseDimension      = 0;                              % process noise dimension
    model.observationNoiseDimension  = 7;                              % observation noise dimension
    model.params                     = initArgs.initialParams;         % setup parameter vector buffer
    model                            = setparams(model, ...
                                            initArgs.initialParams, ...
                                            initArgs.xRaySources, ...
                                            initArgs.earthEphemeris, ...
                                            initArgs.sunEphemeris);

    % Setup process noise source
    processNoiseArg.type           = 'gaussian';
    processNoiseArg.covarianceType = 'full';
    processNoiseArg.tag            = 'GSSM observation noise observation';
    processNoiseArg.dimension      = model.processNoiseDimension;
    processNoiseArg.mean           = 0;
    processNoiseArg.covariance     = 0;
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
function updatedModel = setparams(model, params, xRaySources, earthEphemeris, sunEphemeris)
    % Function to unpack a column vector containing system parameters into specific forms
    % needed by FFUN, HFUN and possibly defined sub-functional objects. Both the vectorized (packed)
    % form of the parameters as well as the unpacked forms are stored within the model data structure.
    updatedModel = deepClone(model);
    updatedModel.params = params;
    
    updatedModel.timeTillCurrentEpoch   = params(1);
    updatedModel.sampleTime             = params(2);
    updatedModel.time                   = params(3);
    updatedModel.xRaySources            = xRaySources;
    updatedModel.earthEphemerisX        = earthEphemeris(1);
    updatedModel.earthEphemerisY        = earthEphemeris(2);
    updatedModel.earthEphemerisZ        = earthEphemeris(3);
    updatedModel.sunEphemerisX          = sunEphemeris(1);
    updatedModel.sunEphemerisY          = sunEphemeris(2);
    updatedModel.sunEphemerisZ          = sunEphemeris(3);
end

function newState = ffun(model, state, noise, stateControl)
    % State transition function (system dynamics).    
%     [~, tmp] = ode113( @(t,y) EquationOfMotion(t, ...
%             y, ...
%             model.acceleration, ...
%             model.angularVelocity, ...
%             model.timeTillCurrentEpoch, ...
%             model.sampleTime ), ...
%         [model.time - model.sampleTime, model.time], ...
%         state ...
%     );
%     newState = tmp(end, :);
    newState = rungeKuttaFourOrderWithFixedStep(...
        @(t,y) equationOfMotionFreeFly(t, ...
            y, ...
            model.timeTillCurrentEpoch), ...
        state, ...
        model.time, ...
        model.sampleTime...
    );
    
    if ~isempty(noise)
        newState = newState + noise';
    end
    
    if ~isempty(stateControl)
        newState = newState + stateControl;
    end
end

function observ = hfun(gssModel, state, noise, observationControl) % first argument is a model.
    % State observation function.
    earthEphemeris.x = gssModel.earthEphemerisX;
    earthEphemeris.y = gssModel.earthEphemerisY;
    earthEphemeris.z = gssModel.earthEphemerisZ;
    
    sunEphemeris.x = gssModel.sunEphemerisX;
    sunEphemeris.y = gssModel.sunEphemerisY;
    sunEphemeris.z = gssModel.sunEphemerisZ;
    
    diffToa = calculateDiffToa(gssModel.xRaySources, earthEphemeris, sunEphemeris, state(1:3)');
    observ = diffToa2phase(gssModel.xRaySources, diffToa);
    
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
%   output of hfun, i.e. observ (the predicted system observation) and an actual 'real world' observation OBS.
    innov = observation - predictedObservation;
end