classdef IntegratedInsXRayNS < BaseIntegratedIns
    % IntegratedInsXRayNS. integrated inertial navigation system (INS) and X-Ray navigation system (navigation system by signal of X-Ray sources).
    % Provides method to solve navigation problem via INS and X-Ray nav system.
    
    properties(Access = private)
        xRayNav;
        xRayState;
        initialXRayState;
        intialXRayCov;
        xRayEstimator;
        xRayTimeData;
    end
    
    methods (Access = public)
        function obj = IntegratedInsXRayNS(ins, xRayNav, timeData, initArgs, xRayTimeData, initialXRayState, intialXRayCov, xRayEstimator, reconciliationTime)
            obj.ins = ins;
            obj.xRayNav             = xRayNav;
            obj.timeData            = timeData;
            obj.initArgs            = initArgs;
            obj.initialXRayState    = initialXRayState;
            obj.intialXRayCov       = intialXRayCov;
            obj.xRayEstimator       = xRayEstimator;
            obj.xRayTimeData        = xRayTimeData;
            obj.xRayState           = [];
            obj.reconciliationTime  = reconciliationTime;
        end
    end
    
    methods (Access = protected)
        function [ state ] = evaluateSecondaryState(this, ~, ~, timeEnd)
            narginchk(4, 4);
            
            endSample = this.timeData.evalSampleFromTime(timeEnd);
            
            if isempty(this.xRayState)
                this.evaluateXRayState();
            end
            
            state = this.xRayState(:, endSample);
        end
        
        function decompCov = updateFilterParams(this, cov, estimatorType)
            decompCov = [];
            switch estimatorType{1}
                case 'ukf'
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                case 'srukf'
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    decompCov = chol(cov, 'lower');
                case 'cdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    this.inferenceModel.spkfParams = sqrt(25); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCov = chol(cov, 'lower');
                case {'sckf', 'fdckf'}
                    decompCov = chol(cov, 'lower');
                case 'pf'
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gspf'
                    this.inferenceModel.estimateType = 'mean';
                case 'sppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gmsppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
                    alpha = 0.75;   % scale factor (UKF parameter) 1e-3
                    beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
                    kappa = 39;     % 0 is optimal for state dimension = 2 (UKF parameter)
                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'cqkf'
                    this.inferenceModel.cqkfParams = 9; % order of laguerre polynomial
                case 'ghqf'
                    this.inferenceModel.ghkfParams = 1; % order of gauss-hermite polynomial
                case 'sghqf'
                    this.inferenceModel.sghkfParams = [3 3]; % order of gauss-hermite polynomial & manner
                otherwise
                    % do nothing by default
            end
        end
        
        function [ description ] = tag(~)
            description = 'State estimation for loosely coupled Ins & X-Ray integrated system';
        end
    end
    
    methods (Access = private)
        function evaluateXRayState(this)
            this.xRayState = this.xRayNav.resolve(this.initialXRayState, this.intialXRayCov, this.xRayEstimator, 0);
        end
    end
    
end
