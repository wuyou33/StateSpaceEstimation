function inferenceDataStructure = inferenceDataGenerator( args )
    % inferenceDataGenerator. Generate inference data structure from a generalized state space model and user defined inference parameters.
    %
    %   INPUT
    %       args.type         Inference (estimation) type : 'state', 'parameter' or 'joint';
    %           .tag          Arbitrary user defined ID tag string;
    %           .model        Generalized state space model descriptor
    %           .updateType   Update type : Time-and-Measurement Update is the default.
    %
    %%
    if (nargin < 1)
        error('[ inferenceDataGenerator ] Not enough inputs.');
    end
    
    if ~isstruct(args)
        error('[ inferenceDataGenerator::args ] Input argument must be an argument data structure args.');
    end
    
    if ~ischar(args.type)
        error('[ inferenceDataGenerator::args ] args.type must be a string indicating the inference type, i.e. ''state'', ''parameter'' or ''joint''.');
    end
    
    %%
    inferenceDataStructure.type             = 'InferenceDataStructure';
    inferenceDataStructure.inferenceType    = args.type;
    inferenceDataStructure.model            = args.model;
    
    if isfield(args, 'tag')
        inferenceDataStructure.tag  = args.tag;
    else
        warning('[inferenceDataGenerator::args] set empty tag by default. please specify tag');
        inferenceDataStructure.tag  = '';
    end
    
    %%
    switch args.type
        case 'state'
            inferenceDataStructure.stateDimension             = args.model.stateDimension;
            inferenceDataStructure.observationDimension       = args.model.observationDimension;
            inferenceDataStructure.controlInputDimension      = args.model.controlInputDimension;     % exogenous input 1 dimension
            inferenceDataStructure.control2InputDimension     = args.model.control2InputDimension;    % exogenous input 2 dimension
            inferenceDataStructure.processNoiseDimension      = args.model.processNoise.dimension;
            inferenceDataStructure.observationNoiseDimension  = args.model.observationNoise.dimension;
            inferenceDataStructure.stateTransitionFun         = @stateTransitionFun;
            inferenceDataStructure.stateObservationFun        = @stateObservationFun;
            inferenceDataStructure.stateTransitionPriorFun    = @stateTransitionPriorFun;
            inferenceDataStructure.likelihoodStateFun         = @likelihoodStateFun;
            inferenceDataStructure.innovationModelFunc        = @innovationModelFunc;
            inferenceDataStructure.linearize                  = @linearize;
            
            % Index vectors indicating the presence of angular components in the state and observation vectors
            if isfield(args.model, 'stateAngleCompIdxVec'),
                inferenceDataStructure.stateAngleCompIdxVec = args.model.stateAngleCompIdxVec;
            else
                inferenceDataStructure.stateAngleCompIdxVec = [];
            end
            
            if isfield(args.model, 'obsAngleCompIdxVec'),
                inferenceDataStructure.obsAngleCompIdxVec = args.model.obsAngleCompIdxVec;
            else
                inferenceDataStructure.obsAngleCompIdxVec = [];
            end
            
        case 'parameter'
            error('[ inferenceDataGenerator::args::type ] parameter estimation not implelemented');
            
        case 'joint'
            error('[inferenceDataGenerator::args::type ] joint estimation not implelemented');
            
        otherwise
            error(['[ inferenceDataGenerator::args::type ] Inference type ''' args.type ''' not supported.']);
    end
    
    inferenceDataStructure.updateType = 'TMU';
end
%%
function newState = stateTransitionFun(inferenceDataModel, state, stateNoise, control)
    %  State transition function of meta system for state estimation
    
    newState = inferenceDataModel.model.stateTransitionFun(inferenceDataModel.model, state, stateNoise, control);
end
%%
function observation = stateObservationFun(inferenceDataModel, state, observationNoise, control)
    %  State observation function of meta system for state estimation
    
    observation = inferenceDataModel.model.stateObservationFun( inferenceDataModel.model, state, observationNoise, control);
end
%%
function tranPrior = stateTransitionPriorFun(inferenceDataModel, predictedState, state, control, processNoiseDataSet)
    %  Calculates the transition prior probability P(x_k|x_(k-1))
    
    tranPrior = inferenceDataModel.model.stateTransitionPriorFun(inferenceDataModel.model, predictedState, state, control, processNoiseDataSet);
end
%%
function llh = likelihoodStateFun(inferenceDataModel, observation, state, control, observationNoiseDS)
    % Calculates the likelood of a real-world observation obs given
    % a realization of the predicted observation for a given state,
    % i.e. p(z|x) = p(observation|state)
    
    llh = inferenceDataModel.model.observationLikelihoodFun(inferenceDataModel.model, observation, state, control, observationNoiseDS);
end
%%
function innov = innovationModelFunc(inferenceDataModel, observation, predictedObservation)
    %  innovationModelFunc. Calculates the innovation signal (difference) between the
    %   output of stateObservationFun, i.e. predictedObservation (the predicted system observation) and an actual
    %   'real world' observation observation. This function might be as simple as
    %   innov = observation - predictedObservation, which is the default case, but can also be more
    %   complex for complex measurement processes where for example multiple (possibly false)
    %   observations can be observed for a given hidden ground truth.
    
    innov = inferenceDataModel.model.innovationModelFunc(inferenceDataModel.model, observation, predictedObservation);
end
%%
function varargout = linearize(inferenceDataModel, state, stateNoise, observNoise, control1, control2, varargin)
    %  linearize. Linearization function of meta system for state estimation.
    %
    %    out = linearize(inferenceDataModel, state, stateNoise, observNoise, control1, control2, varargin)
    %
    %    INPUT
    %         inferenceDataModel    inference data structure;
    %         state                 meta system state vector;
    %         stateNoise            meta system process noise vector;
    %         observNoise           meta system observation noise vector;
    %         control1              meta system exogenous input 1;
    %         control2              meta system exogenous input 2;
    %         varargin              linearization terms wanted, e.g. 'A', 'B', 'G', ... .
    %
    %    OUTPUT
    %         varargout     linearization terms corresponding with varargin strings.
    
    numRequestedTerm = length(varargin);
    varargout = cell(numRequestedTerm, 1);
    
    for k = 1:numRequestedTerm,
        varargout{k} = inferenceDataModel.model.linearize(inferenceDataModel.model, state, stateNoise, observNoise, control1, control2, varargin{k});
    end
end
%%