classdef SnsImitator < handle
    %SNSIMITATOR: Imitator of satellity navigation system. If real sns
    % receiber is available, then the receiber will used, otherwise
    % real state with additive white noise will used.   
    
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
        function obj = SnsImitator(state, noiseType) 
            if nargin >= 1
                num = length(state);
                if nargin == 2 
                    switch noiseType
                        case 'gaussian'
                            noiseTraj = 2e-2*randn(3, num); % exprnd(ones(3, dim)) % km
                            noiseVel  = 1e-5*randn(3, num); % exprnd(ones(3, dim)) % km / sec
                        case 'exp'
                            noiseTraj = 1e-2*exprnd(ones(3, num));
                            noiseVel  = 1e-5*exprnd(ones(3, num));
                        case 'wiener'
                            noiseProc = WienerProcess([0; 0; 0], [1; 1; 1]);
                            noiseTraj = 1e-2*noiseProc.simulate(1e-3, num);
                            noiseVel  = 1e-5*noiseProc.simulate(1e-3, num);
                        case 'markov1'
                            noiseProc = ExpCorrGaussianProcess(10*ones(3, 1), zeros(3, 1), ones(3, 1));
                            noiseTraj = 1e-2*noiseProc.simulate(num);
                            noiseVel  = 1e-5*noiseProc.simulate(num);
                        otherwise
                            error('[ SnsImitator ] unknown noiseType');
                    end
                else
                    noiseTraj = 1e-2*randn(3, num); % exprnd(ones(3, dim))
                    noiseVel  = 1e-5*randn(3, num); % exprnd(ones(3, dim))
                end
                
                obj.trajectory = state(1:3, :) + noiseTraj;
                obj.velocity   = state(4:6, :) + noiseVel;
                obj.trueTrajectory = state(1:3, :);
                obj.trueVelocity   = state(4:6, :);
            else
                gnssStateEstimated = load('exclude/GnssState.mat');            
                obj.trajectory = gnssStateEstimated.GNSS_state2(1:3, :);
                obj.velocity   = gnssStateEstimated.GNSS_state2(4:6, :);
                
                gnssStateTrue = load('exclude/ScEphem.mat');
                obj.trueTrajectory = gnssStateTrue.SC2(1:3, :);
                obj.trueVelocity   = gnssStateTrue.SC2(4:6, :);
            end
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
