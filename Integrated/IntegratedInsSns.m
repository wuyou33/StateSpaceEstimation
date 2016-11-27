classdef IntegratedInsSns < BaseIntegratedIns
    % IntegratedInsSns. integrated inertial navigation system (INS) and satellity navigation system (SNS).
    % Provides method to solve navigation problem via INS and SNS (GPS)
    
    properties(Access = private)
        sns;
        snsTimeData;
    end
    
    methods (Access = public)
        function obj = IntegratedInsSns(ins, sns, timeData, initArgs, snsTimeData, reconciliationTime)
            obj.ins                 = ins;
            obj.sns                 = sns;
            obj.timeData            = timeData;
            obj.initArgs            = initArgs;
            obj.snsTimeData         = snsTimeData;
            obj.reconciliationTime  = reconciliationTime;
        end
    end
    
    methods (Access = protected)
        function [ state ] = evaluateSecondaryState(this, ~, ~, timeEnd)
            endSample = this.snsTimeData.evalSampleFromTime(timeEnd);
            state = this.sns.getState(endSample);
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
                    decompCov = svdDecomposition(cov);
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
            description = 'State estimation for loosely coupled Ins & Sns integrated system';
        end
    end
    
end
