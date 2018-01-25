classdef IntegratedIns < BaseIntegratedIns
    % IntegratedIns. Integrated inertial navigation system (INS) and abstract navigation system which provide information about position and velocity.
    % Provides method to solve navigation problem via INS and any other nav system which provide information about position and velocity.
    
    properties(Access = private)
        fns; % abstract navigation system which can provide information about position and velocity.
        fnsTimeData;
    end
    
    methods (Access = public)
        function obj = IntegratedIns(ins, fns, timeData, initArgs, fnsTimeData, reconciliationTime)
            obj.ins                 = ins;
            obj.fns                 = fns;
            obj.timeData            = timeData;
            obj.initArgs            = initArgs;
            obj.fnsTimeData         = fnsTimeData;
            obj.reconciliationTime  = reconciliationTime;
        end
    end
    
    methods (Access = protected)
        function [ state ] = evaluateSecondaryState(this, ~, ~, timeEnd)
            endSample = this.fnsTimeData.evalSampleFromTime(timeEnd);
            state = this.fns.getState(endSample);
        end
        
        function decompCov = updateFilterParams(this, cov, estimatorType)
            decompCov = [];
            alpha = 0.75;   % scale factor (UKF parameter) 1e-3
            beta  = 0.85;   % 2 is a optimal setting for Gaussian priors (UKF parameter)
            kappa = 22;     % 0 is optimal for state dimension = 2 (UKF parameter)
            
            switch estimatorType{1}
                case 'ukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                case 'srukf'
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    decompCov = chol(cov, 'lower');
                case 'cdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                case 'srcdkf'
                    this.inferenceModel.spkfParams = sqrt(7); % scale factor (CDKF parameter h) default sqrt(3)
                    decompCov = chol(cov, 'lower');
                case {'sckf', 'fdckf'}
                    decompCov = chol(cov', 'lower');
                case 'pf'
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gspf'
                    this.inferenceModel.estimateType = 'mean';
                case 'sppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');                    
                    this.inferenceModel.spkfParams = [alpha beta kappa];
                    this.inferenceModel.resampleThreshold   = 1;
                    this.inferenceModel.estimateType        = 'mean';
                case 'gmsppf'
                    this.inferenceModel.spkfType    = 'srukf';
                    decompCov                       = chol(cov, 'lower');
                    
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
            description = 'State estimation for loosely coupled Ins & ans integrated system';
        end
    end
    
end
