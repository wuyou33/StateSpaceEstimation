classdef AccelerationInBodyFrame < handle
    % Provide acceleration cuased by external forces.
    % [km / sec^2]
    %    
    properties (Access = private)
        acceleration; % [km / sec^2]
        model;
        simulationNumber;
        dt
    end
    
    properties (Dependent)
        % Acceleration in body frame [km / sec^2].
        Acceleration;
    end
    
    properties (Constant, Access = private)
        Dimension = 3;
    end
    
    methods        
        function obj = AccelerationInBodyFrame(simulationNumber, dt, mu, sigma)
            %% Constructor. 
            %   mu               [km / sec^2]
            %   sigma            [km / sec^2]
            %   simulationNumber [-]
            %   dt               [sec]
            
            % bmModel = bm(mu, sigma);
            % bmModel.StartState = 0;
            % obj.acceleration = simulate(bmModel, simulationNumber);
            % obj.model = bmModel;
            obj.simulationNumber = simulationNumber;
            if (nargin == 1)
                obj.dt    = NaN;
                obj.model = NaN;
                obj.acceleration = zeros(simulationNumber, AccelerationInBodyFrame.Dimension);
            else
                obj.dt = dt;                
                dw = WienerProcess(mu, sigma);
                obj.model = dw;
                obj.acceleration = dw.simulate(dt, simulationNumber);
            end
        end
        
        function val = get.Acceleration(this)
            val = this.acceleration;
        end
    end
    
    methods(Static, Access = public)
        function obj = Empty(simulationNumber)
            obj = AccelerationInBodyFrame(simulationNumber);
        end
    end
end