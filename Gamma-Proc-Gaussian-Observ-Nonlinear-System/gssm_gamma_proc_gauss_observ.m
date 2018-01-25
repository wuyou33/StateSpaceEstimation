function [ varargout ] = gssm_gamma_proc_gauss_observ( func, varargin )
    % GSSM_GAMMA_PROC_GAUSS_OBSERV
    % The model is a scalar nonlinear system with Gamma process and Gaussian observation noise.
    %
    switch func
        case 'init'
            varargout{1} = init(varargin);
        otherwise
            error(['Function ''' func ''' not supported.']);
    end
end

function model = init(init_args)
    model.type                       = 'gssm';
    model.tag                        = 'nonlinear system with Gamma process and Gaussian observation noise';
    
    model.stateTransitionFun         = @ffun;        % function handle to state transition function
    model.stateObservationFun        = @hfun;        % function handle to state observation function
    model.stateTransitionPriorFun    = @prior;       % function handle to the state transition function that calculates P(x(k)|x(k-1)),
    model.observationLikelihoodFun   = @likelihood;  % function handle to the observation likelihood function that calculates p(y(k)|x(k)),
    model.innovationModelFunc        = @innovation;  % Function-handle to the innovation model function that calculates the difference between the output
    % of the observation function (hfun) and the actual 'real-world' measurement/observation of that signal
    model.linearize                  = @linearize;                     % Function-handle to the linearization function that calculates Jacobians e.t.c.
    model.setParams                  = @set_params;   % function handle to SETPARAMS
    
    model.stateDimension             = 1;
    model.observationDimension       = 1;
    model.paramDimension             = 2;
    model.controlInputDimension      = 1;
    model.control2InputDimension     = 1;
    model.processNoiseDimension      = 1;
    model.observationNoiseDimension  = 1;
    
    % Setup process noise source
    processNoiseArg.type        = 'gamma';
    processNoiseArg.tag         = 'GSSM process (state) noise';
    processNoiseArg.dimension   = model.processNoiseDimension;
    processNoiseArg.alpha       = 3;
    processNoiseArg.beta        = 0.5;
    
    model.processNoise = generateNoiseDataSet(processNoiseArg);
    
    
    % Setup observation noise source
    observationNoiseArg.type           = 'gaussian';
    observationNoiseArg.covarianceType = 'full';
    observationNoiseArg.tag            = 'GSSM observation noise';
    observationNoiseArg.dimension      = model.observationNoiseDimension;
    observationNoiseArg.mean           = 0;
    observationNoiseArg.covariance     = 1e-5;
    
    model.observationNoise = generateNoiseDataSet(observationNoiseArg);
    
    model.params = zeros(model.paramDimension, 1);
    model = set_params(model, [4e-2 0.5]); % [omega phi]
end

function model = set_params(model, params, index_vector)
    switch(nargin)
        case 2
            model.params = params(:);
        case 3
            model.params(index_vector) = params(:);
        otherwise
            error('[ set_params ] Incorrect number of input arguments.');
    end
end

function new_state = ffun(model, state, v, u1)
    new_state = 1 + sin(model.params(1) *pi .* u1) + model.params(2) * state;
    
    if (~isempty(v))
        new_state = new_state + v;
    end
end

function observ = hfun(model, state, obs_noise, u2)
    [~, obs_count] = size(state);
    observ = zeros(model.observationDimension, obs_count);
    
    for i = 1 : obs_count
        observ(i) = model.params(2) * state(:, i) .^ 2;
    end
    
    if (~isempty(obs_noise))
        observ = observ + obs_noise;
    end
end

function tranprior = prior(model, next_state, state, u1, proc_noise_infer_model)
    proc_deviation = next_state - ffun(model, state, [], u1);
    tranprior = proc_noise_infer_model.likelihood(proc_noise_infer_model, proc_deviation);
end

function llh = likelihood(model, obs, state, u2, obs_noise_infer_model)
    obs_deviation = obs - hfun(model, state, [], u2);
    llh = obs_noise_infer_model.likelihood(obs_noise_infer_model, obs_deviation);
end

function innov = innovation(~, real_obs, predicted_observ)
    innov = real_obs - predicted_observ;
end

function out = linearize(model, state, proc_noise, obs_noise, u1, u2, term, index_vector)
    if (nargin < 7)
        error('[ linearize ] Not enough input arguments!');
    end
    
    switch (term)
        case 'A'
            % A = df / dstate
            out = model.params(2);
        case 'B'
            % B = df / du1
            out = [];
        case 'C'
            % C = dh / dx
            if (u2 <= 30)
                out = 2 * model.params(2) * state;
            else
                out = model.params(2);
            end
        case 'D'
            % D = dh / du2
            out = [];
        case 'G'
            %G = df / dv
            out = 1;
        case 'H'
            % H = dh / dn
            out = 1;
        case 'JFW'
            % dffun / dparameters
            out = [ cos(model.params(1) * pi * u1) * pi * u1 state];
        case 'JHW'
            % dhfun / dparameters
            if (u2 <= 30)
                out = [0 state^2];
            else
                out = [0 state];
            end
        otherwise
            error('[ linearize ] Invalid model term requested!');
    end
    
    if (nargin == 8)
        out = out(:, index_vector);
    end
end