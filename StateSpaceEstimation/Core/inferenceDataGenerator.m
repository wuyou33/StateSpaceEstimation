function inferenceDataStructure = inferenceDataGenerator( args )
% Generate inference data structure from a generalized state space model and user defined inference parameters.
% INPUT ARGUMENT DATA STRUCTURE:
% args.type               : (string)      Inference (estimation) type : 'state', 'parameter' or 'joint'
%     .tag                : (string)      Arbitrary user defined ID tag string
%     .model              : (gssm)        Generalized state space model descriptor
%     .updateType         : (string)      Update type : Time-and-Measurement Update is the default. Other options are : TU - time update only or MU - measurement update only.
%%
        if (nargin < 1)
            error(' [ inferenceDataGenerator ] Not enough inputs.');
        end

        if ~isstruct(args)
            error(' [ inferenceDataGenerator ] Input argument must be an argument data structure ArgDS.');
        end

        if ~ischar(args.type)
            error(' [ inferenceDataGenerator ] ArgDS.type must be a string indicating the inference type, i.e. ''state'', ''parameter'' or ''joint''.');
        end

    %%
        inferenceDataStructure.type             = 'InferenceDataStructure';
        inferenceDataStructure.inferenceType    = args.type;
        inferenceDataStructure.model            = args.model;

        if isfield(args,'tag'),
            inferenceDataStructure.tag  = args.tag;
        else
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
            inferenceDataStructure.innovationState            = @innovationState;

            % Index vectors indicating the presence of angular components in the state and observation vectors
            if isfield(args.model,'stateAngleCompIdxVec'),
                inferenceDataStructure.stateAngleCompIdxVec = args.model.stateAngleCompIdxVec;
            else
                inferenceDataStructure.stateAngleCompIdxVec = [];
            end

            if isfield(args.model,'obsAngleCompIdxVec'),
                inferenceDataStructure.obsAngleCompIdxVec = args.model.obsAngleCompIdxVec;
            else
                inferenceDataStructure.obsAngleCompIdxVec = [];
            end

        case 'parameter'
            error(' [ inferenceDataGenerator ] parameter estimation not implelemented');

        case 'joint'
            error(' [inferenceDataGenerator ] joint estimation not implelemented');

        otherwise
            error([' [ inferenceDataGenerator ] Inference type ''' args.type ''' not supported.']);        
    end

    inferenceDataStructure.updateType = 'TMU';
end

function newState = stateTransitionFun(inferenceDS, state, stateNoise, control)
    %  State transition function of meta system for state estimation

    newState = inferenceDS.model.stateTransitionFun( inferenceDS.model, state, stateNoise, control);
end

function observation = stateObservationFun(inferenceDS, state, observationNoise, control)
    %  State observation function of meta system for state estimation
    
    observation = inferenceDS.model.stateObservationFun( inferenceDS.model, state, observationNoise, control);
end

function tran_prior = stateTransitionPriorFun(inferenceDS, predictedState, state, control, processNoiseDataSet)
    %  Calculates the transition prior probability P(x_k|x_(k-1))
 
    tran_prior = inferenceDS.model.stateTransitionPriorFun( inferenceDS.model, predictedState, state, control, processNoiseDataSet);
end

function llh = likelihoodStateFun(inferenceDS, observation, state, control, observationNoiseDS)
    % Calculates the likelood of a real-world observation obs given
    % a realization of the predicted observation for a given state,
    % i.e. p(y|x) = p(obs|state)
    
    llh = inferenceDS.model.observationLikelihoodFun( inferenceDS.model, observation, state, control, observationNoiseDS);
end

function innov = innovationState(inferenceDS, observation, predictedObservation)
    %  INNOVATION_STATE  Calculates the innovation signal (difference) between the
    %   output of HFUN, i.e. predictedObservation (the predicted system observation) and an actual
    %   'real world' observation observation. This function might be as simple as
    %   INNOV = observation - predictedObservation, which is the default case, but can also be more
    %   complex for complex measurement processes where for example multiple (possibly false)
    %   observations can be observed for a given hidden ground truth.

    innov = inferenceDS.model.innovationModelFunc( inferenceDS.model, observation, predictedObservation);
end