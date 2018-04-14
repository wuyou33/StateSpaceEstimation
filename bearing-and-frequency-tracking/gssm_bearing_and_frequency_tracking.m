function [ varargout ] = gssm_bearing_and_frequency_tracking( func, varargin )
    % General state space model for Bearings and Frequency Tracking of a randomly maneuvering
    %           target relative to a moving observer (submarine).
    %
    %   The following state space model is used :
    %
    %     X(k) = |1 4 0 0 0| X(k-1) + | 8  0  0| V(k-1)
    %            |0 1 0 0 0|          | 4  0  0|
    %            |0 0 1 4 0|          | 0  8  0|
    %            |0 0 0 1 0|          | 0  4  0|
    %            |0 0 0 0 1|          | 0  0  1|
    %    (the measures are given by the sonar every 4 seconds)
    %
    %   Where the state vector is defined as
    %            - the 2D position and velocity vector of the target (relative to a fixed external reference frame),
    %            - the pure tone frequency emitted by the target (very stable).
    %
    %     X(k) = |x1(k)| = |      x-position at time k     |
    %            |x2(k)|   |      x-velocity at time k     |
    %            |x3(k)|   |      y-position at time k     |
    %            |x4(k)|   |      y-velocity at time k     |
    %            |x5(k)|   | pure tone frequency at time k |
    %
    %   And the observations at time k, O(k) are :
    %           - the bearing angle (in radians) from the moving observer (submarine) towards the target,
    %           - the doppler-shifted frequency tone tracked by the observer (submarine).
    %
    %     O(k) = |                bearing = arctan((x3(k)-sub3(k))/(x1(k)-sub1(k)))                            |  +  | v1(k) |
    %            |   frequency = x5(k)*(1+1/1500*((x2(k)-sub2(k))*cos(bearing)+(x4(k)-sub4(k))*sin(bearing)))  |     | v2(k) |
    %
    %   c = 1500 m/s (sound speed in water)
    %
    %   The submarine state is known at each time precisely and described by the following vector :
    %     sub(k) = |sub1(k)| = |  x-position at time k     |
    %              |sub2(k)|   |  x-velocity at time k     |
    %              |sub3(k)|   |  y-position at time k     |
    %              |sub4(k)|   |  y-velocity at time k     |
    %              |sub5(k)|   |  frequency tone at time k | (not used here)
    %
    %   The state dynamics are driven by a 2 dimensional white Gaussian noise source and the
    %   observations are corrupted by additive scalar white Gaussian noise.
    %
    %   See :  Gordon, Salmond & Ewing, "Bayesian State Estimation for Tracking and Guidance Using
    %   the Bootstrap Filter", Journal of Guidance, Control and Dynamics, 1995.
    
    switch func
        case 'init'
            varargout{1} = init(varargin);
        otherwise
            error(['Function ''' func ''' not supported.']);
    end
end

function model = init(init_args)
    
    model.type = 'gssm';
    model.tag  = 'gssm: bearings and frequency tracking';
    
    model.stateDimension              = 5;
    model.observationDimension        = 2;
    model.paramDimension              = 12;
    model.controlInputDimension       = 0;
    model.control2InputDimension      = 5;
    model.processNoiseDimension       = 3;
    model.observationNoiseDimension   = 2;
    
    model.transition_fun    = @ffun;
    model.observation_fun   = @hfun;
    model.prior             = @prior;
    model.likelihood        = @likelihood;
    model.innovation        = @innovation;
    model.set_params        = @set_params;
    
    % indicate that the first (and only component) of the observation vector is an angle measured in radians.
    % This is needed so that the SPKF based algorithms can correctly deal with the angular discontinuity at +- pi radians.
    model.obsAngleCompIdxVec = 1;
    
    
    processNoiseArg.type = 'gaussian';
    processNoiseArg.dimension = model.processNoiseDimension;
    processNoiseArg.mean  = zeros(processNoiseArg.dimension, 1);
    processNoiseArg.covarianceType = 'full';
    processNoiseArg.covariance   = [((1e-3)^2)*eye(processNoiseArg.dimension - 1) zeros(processNoiseArg.dimension - 1, 1); ...
        zeros(1, processNoiseArg.dimension - 1) (1e-4)^2];
    model.processNoise = generate_noise_model(processNoiseArg);
    
    observationNoiseArg.type = 'gaussian';
    observationNoiseArg.dimension = model.observationNoiseDimension;
    observationNoiseArg.mean = zeros(observationNoiseArg.dimension, 1);
    observationNoiseArg.covarianceType ='full';
    observationNoiseArg.covariance  = [0.0175^2 0; 0 0.06^2];
    model.observationNoise = generate_noise_model(observationNoiseArg);
    
    model.params = zeros(model.paramDimension, 1);
    model.A = zeros(model.stateDimension, model.stateDimension);
    model.G = zeros(model.stateDimension, model.processNoiseDimension);
    
    model = set_params(model, [1 4 1 1 4 1 1 8 4 8 4 1]');
end

function model = set_params(model, params, index_vector)
    if (nargin==2)
        model.params = params(:);
    elseif (nargin == 3)
        model.params(index_vector) = params(:);
    else
        error('[ set_params ] Incorrect number of input arguments.');
    end
    
    model.A([1 6 7 13 18 19 25]) = params(1:7);
    model.G([1 2 8 9 15]) = params(8:12);
    
    G = model.G;
    % conversion matrix needed to calculate state transition prior
    model.convFact1 = (G'*G) \ (G');
end

function new_state = ffun(model, state, v, u1)
    if isempty(v)
        new_state = model.A*state;
    else
        new_state = model.A*state + model.G*v;
    end
end

function observ = hfun(model, state, obs_noise, u2)
    observ_ = zeros(2, size(state, 2));
    observ_(1, :) = atan2(state(3, :) - u2(3, :), state(1, :) - u2(1, :));
    observ_(2, :) = state(5, :) .* (1 + 1/1500 * ((u2(2, :) - state(2, :)) .* cos(observ_(1,:))+(u2(4,:)-state(4,:)).*sin(observ_(1,:))));
    
    % Now add the measurement noise... taking care with the discontinueties at +-pi radians
    if isempty(obs_noise),
        observ = observ_;
    else
        observ = observ_ + obs_noise;
        observ(1,:) = fix_angular_discontinuety(observ_(1,:), obs_noise(1,:));
    end
end

function tranprior = prior(model, nextstate, state, u1, proc_noise_model)
    v = model.convFact1 * (nextstate - model.A*state);
    tranprior = proc_noise_model.likelihood( proc_noise_model, v);
end

function llh = likelihood(model, obs, state, u2, observ_noise_model)
    observ =  hfun(model, state, [], u2);
    obs_deviation = obs - observ;
    obs_deviation(1, :) = sub_angle(obs(1,:), observ(1, :));
    
    llh = observ_noise_model.likelihood(observ_noise_model, obs_deviation);
end

function innov = innovation(model, obs, observ)
    innov = obs - observ;
    
    % deal with the discontinueties at +-pi radians
    innov(1, :) = sub_angle(obs(1, :), observ(1, :));
end
