% Generalized state space model for X-Ray navigation system
%%
function [varargout] = gssm_x_nav(func, varargin)
    switch func
        case 'init'
            varargout{1} = init(varargin{1});
        otherwise
            error(['Function ''' func ''' not supported.']);
    end
end
%%
function model = init(args)
    model.type              = 'gssm';
    model.tag               = 'X-Ray Navigation System';
    
    model.set_params        = @set_params; % function handle to set parameters.
    model.transition_fun    = @f_fun; % function handle to state transition function
    model.observation_fun   = @h_fun; % function handle to state observation function
    model.linearize         = @linearize; % Function-handle to the linearization function that calculates Jacobians e.t.c.
    model.prior             = @prior; % function handle to the state transition function that calculates P(x(k)|x(k-1)), is a requirement for particle filter
    model.likelihood        = @likelihood; % function handle to the observation likelihood function that calculates p(y(k)|x(k)), is a requirement for particle filter
    model.innovation        = @innovation; % Function-handle to the innovation model function that calculates the difference between the output
    % of the observation function (hfun) and the actual 'real-world' measurement/observation of that signal,
    % is a requirement for particle filter
    
    model.stateDimension             = 6;
    model.observationDimension       = 4;
    model.paramDimension             = 0; % param
    model.controlInputDimension      = 0; % exogenous control input 1 dimension
    model.control2InputDimension     = 0; % exogenous control input 2 dimension
    model.processNoiseDimension      = 6; % process noise dimension
    model.observationNoiseDimension  = 4; % observation noise dimension
    
    model = set_params(model, args.current_epoch, args.current_time);
    model = set_permanent_params(model, args.mass, args.gravity_model, args.sample_time, args.photon_detector, args.earth_ephemeris, args.sun_ephemeris);
    
    % Setup process noise source
    process_noise_arg.type           = 'gaussian';
    process_noise_arg.covarianceType = 'full';
    process_noise_arg.tag            = 'GSSM state noise';
    process_noise_arg.dimension      = model.processNoiseDimension;
    process_noise_arg.mean           = zeros(model.processNoiseDimension, 1);
    process_noise_arg.covariance     = [(1e-5*eye(3)).^2 zeros(3); zeros(3) (1e-9*eye(3)).^2]; % [km^2  n/a; n/a (km/sec)^2]
    
    model.processNoise = generate_noise_model(process_noise_arg);
    
    
    % Setup observation noise source
    observation_noise_arg.type           = 'gaussian';
    observation_noise_arg.covarianceType = 'full';
    observation_noise_arg.tag            = 'GSSM observation noise';
    observation_noise_arg.dimension      = model.observationNoiseDimension;
    observation_noise_arg.mean           = zeros(model.observationNoiseDimension, 1);
    observation_noise_arg.covariance     = model.photon_detector.dtoa_cov;
    
    model.observationNoise = generate_noise_model(observation_noise_arg);
end
%%
function updated_model = set_permanent_params(model, mass, gravity_model, sample_time, photon_detector, earth_ephemeris, sun_ephemeris)
    % Initializing permanent properties of model which are constant, at least, during one run.
    updated_model = deep_clone(model);
    
    updated_model.mass              = mass;
    updated_model.gravity_model     = gravity_model;
    updated_model.sample_time       = sample_time;
    updated_model.photon_detector   = photon_detector;
    updated_model.earth_ephemeris   = earth_ephemeris;
    updated_model.sun_ephemeris     = sun_ephemeris;
end

function updated_model = set_params(model, current_epoch, current_time)
    % Updating transient properties of model which may be changed at time to time.
    updated_model = deep_clone(model);
    
    updated_model.current_epoch     = current_epoch;
    updated_model.current_time      = current_time;
end

function new_state = f_fun(model, state, noise, state_control)
    % state transition function (system dynamics).
    t_span = [model.current_time - model.sample_time; model.current_time];
    ode_fun = @(t, y) equation_of_motion_free_fly(t, y, model.current_epoch, model.gravity_model, model.mass);
    
    new_state = zeros(size(state));
    count = size(state, 2);
    for i = 1:count
        [~, tmp]        = ode_euler(ode_fun, t_span, state(:, i), model.sample_time);
        new_state(:, i) = tmp(end, :)';
    end
    
    if ~isempty(noise)
        new_state  = new_state + noise;
    end
    
    if ~isempty(state_control)
        new_state = new_state - column_vector_replicate(state_control, count);
    end
end

function observ = h_fun(model, state, noise, observation_control)
    % Map state to observation function.
    count = size(state, 2);
    
    curr_earth_ephemeris = model.earth_ephemeris(model.current_time);
    earth_ephemeris.x = row_vector_replicate(curr_earth_ephemeris.x, count);
    earth_ephemeris.y = row_vector_replicate(curr_earth_ephemeris.y, count);
    earth_ephemeris.z = row_vector_replicate(curr_earth_ephemeris.z, count);
    
    curr_sun_ephemeris = model.sun_ephemeris(model.current_time);
    sun_ephemeris.x = row_vector_replicate(curr_sun_ephemeris.x, count);
    sun_ephemeris.y = row_vector_replicate(curr_sun_ephemeris.y, count);
    sun_ephemeris.z = row_vector_replicate(curr_sun_ephemeris.z, count);
    
    observ = model.photon_detector.eval_dtoa(state(1:3, :), earth_ephemeris, sun_ephemeris);
    
    if ~isempty(noise)
        observ = observ + noise;
    end
    
    if ~isempty(observation_control)
        observ = observ + observation_control;
    end
end

function trans_prior = prior(model, predicted_state, state, state_control, proc_noise_infer_model)
    %   Calculates P(nextstate|state).
    x = predicted_state - f_fun(model, state, [], state_control);
    trans_prior = proc_noise_infer_model.likelihood(proc_noise_infer_model, x);
end

function llh = likelihood(model, observation, state, observation_control, observ_noise_infer_model)
    %   Observation likelihood function
    z = observation - h_fun(model, state, [], observation_control);
    llh = observ_noise_infer_model.likelihood( observ_noise_infer_model, z);
end

function innov = innovation(~, observation, predicted_observation) % first argument is a model.
    %   Calculates the innovation signal (difference) between the
    %   output of hfun, i.e. observ (the predicted system observation) and an actual 'real world' observation.
    innov = observation - predicted_observation;
end

function out = linearize(model, state, ~, ~, ~, ~, term, ~)
    % (model, state, stateNoise, observNoise, control1, control2, term, index_vector)
    narginchk(7, 8);
    
    switch (term)
        case 'F'
            % A = df / dstate
            func = @(y) equation_of_motion_free_fly(model.time, y, model.timeTillCurrentEpoch, model.gravityModel, model.mass, model.sampleTime, model.startTime);
            [ out, ~ ] = jacobianest(func, state);
        case 'B'
            % B = df / du1, where u1 - control1
            out = zeros(model.stateDimension);
        case 'C'
            % C = dh / dx
            earthEphemeris.x = model.earthEphemerisX;
            earthEphemeris.y = model.earthEphemerisY;
            earthEphemeris.z = model.earthEphemerisZ;
            
            sunEphemeris.x = model.sunEphemerisX;
            sunEphemeris.y = model.sunEphemerisY;
            sunEphemeris.z = model.sunEphemerisZ;
            func = @(y) calculateDiffToa(model.xRaySources, earthEphemeris, sunEphemeris, y);
            [ out, ~ ] = jacobianest(func, state);
        case 'D'
            % D = dh / du2, where u2 - control2
            out = zeros(model.observationDimension);
        case 'G'
            % G = df / dv
            out = eye(model.stateDimension);
        case 'H'
            % H = dh / dn
            out = eye(model.observationDimension);
        case 'JFW'
            % dffun / dparameters
            out = [];
        case 'JHW'
            % dhfun/dparameters
            out = [];
        otherwise
            error('[ linearize::term ] Invalid model term requested!');
    end
end
