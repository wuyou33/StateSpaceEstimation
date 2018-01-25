classdef NavImitator < handle
    % NavImitator. Imitator of navigation system which provide information about velocity and position.
    % Use real state with additive noise (noise model dependent from arguments) will used.
    % Support following noise type:
    %   a) gaussian;
    %   b) exp;
    %   c) wiener;
    %   d) markov1 (1-st order markov process).
    % 
    properties (Access = private)
        trajectory;         % [km]
        velocity;           % [km / sec]
        trueTrajectory;     % [km]
        trueVelocity;       % [km / sec]
    end
    
    properties (Dependent)
        Trajectory;     % [km]
        Velocity;       % [km / sec]
        TrueTrajectory; % [km]
        TrueVelocity;   % [km / sec]
    end
    
    methods (Access = public)
        function obj = NavImitator(state, sigmaTrajectory, meanTrajectory, sigmaVelocity, meanVelocity, noiseType)
            num = length(state);
            
            switch noiseType
                case 'gaussian'
                    noiseTraj = sqrt(sigmaTrajectory / 3)*randn(3, num) + meanTrajectory / 3 * ones(3, num); % km
                    noiseVel  = sqrt(sigmaVelocity / 3)*randn(3, num) + meanVelocity / 3 * ones(3, num); % km / sec
                case 'exp'
                    noiseTraj = sqrt(sigmaTrajectory / 3)*exprnd(ones(3, num)) + meanTrajectory / 3 * ones(3, num);
                    noiseVel  = sqrt(sigmaVelocity / 3)*exprnd(ones(3, num)) + meanVelocity / 3 * ones(3, num);
                case 'wiener'
                    noiseProc = WienerProcess([0; 0; 0], [1; 1; 1]);
                    noiseTraj = sqrt(sigmaTrajectory / 3)*noiseProc.simulate(1e-3, num) + meanTrajectory / 3 * ones(3, num);
                    noiseVel  = sqrt(sigmaVelocity / 3)*noiseProc.simulate(1e-3, num) + meanVelocity / 3 * ones(3, num);
                case 'markov1'
                    noiseProc = ExpCorrGaussianProcess(10*ones(3, 1), zeros(3, 1), ones(3, 1));
                    noiseTraj = sqrt(sigmaTrajectory / 3)*noiseProc.simulate(num) + meanTrajectory / 3 * ones(3, num);
                    noiseVel  = sqrt(sigmaVelocity / 3)*noiseProc.simulate(num) + meanVelocity / 3 * ones(3, num);
                case 'realImitator'
                    error('[ NavImitator ] not supported');
                otherwise
                    error('[ NavImitator ] unknown noiseType');
            end
            
            obj.trajectory = state(1:3, :) + noiseTraj;
            obj.velocity   = state(4:6, :) + noiseVel;
            obj.trueTrajectory = state(1:3, :);
            obj.trueVelocity   = state(4:6, :);
        end
        
        function res = getTrajectory(this, index)
            res = this.trajectory(:, index);
        end
        
        function res = getVelocity(this, index)
            res = this.velocity(:, index);
        end
        
        function res = getState(this, index)
            res = [this.getTrajectory(index); this.getVelocity(index)];
        end
        
        function res = getTrajectoryError(this)
            res = 1e3*(this.trajectory - this.trueTrajectory);
        end
        
        function res = getVelocityError(this)
            res = 1e3*(this.velocity - this.trueVelocity);
        end
    end
    
    methods
        function res = get.Trajectory(this)
            res = this.trajectory;
        end
        
        function res = get.Velocity(this)
            res = this.velocity;
        end
        
        function res = get.TrueTrajectory(this)
            res = this.trueTrajectory;
        end
        
        function res = get.TrueVelocity(this)
            res = this.trueVelocity;
        end
    end
    
end
