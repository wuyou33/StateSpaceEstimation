classdef AccelerationInBodyFrame < handle
    properties (Access = private)
        acceleration;
        model;
        simulationNumber;
        dt
    end
    
    properties (Dependent)
        Acceleration;
    end
    
    properties (Constant, Access=private)
        Dimension = 3;
    end
    
    methods
        function obj = AccelerationInBodyFrame(simulationNumber, dt, mu, sigma)
            %             bmModel = bm(mu, sigma);
            %             bmModel.StartState = 0;
            %             obj.acceleration = simulate(bmModel, simulationNumber);
            %             obj.model = bmModel;
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
    methods(Static, Access=public)
        function obj = Empty(simulationNumber)
            obj = AccelerationInBodyFrame(simulationNumber);
        end
    end
end