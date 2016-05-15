classdef AngularVelocityInBodyFrame < handle
    properties (Access = private)
        angularVelocity;
        model;
        simulationNumber;
        dt
    end
    
    properties (Dependent)
        Velocity;
    end
    
    properties (Constant, Access=private)
        Dimension = 3;
    end
    
    methods
        function obj = AngularVelocityInBodyFrame(simulationNumber, dt, mu, sigma)
            %             bmModel = bm(mu, sigma);
            %             bmModel.StartState = 0;
            %
            %             this.angularVelocity = simulate(bmModel, simulationNumber);
            %             this.model = bmModel;
            obj.simulationNumber = simulationNumber;
            
            if (nargin == 1)
                obj.dt    = NaN;
                obj.model = NaN;
                obj.angularVelocity = zeros(simulationNumber, AngularVelocityInBodyFrame.Dimension);
            else
                obj.dt = dt;                
                dw = WienerProcess(mu, sigma);
                obj.model = dw;
                obj.angularVelocity = dw.simulate(dt, simulationNumber);
            end
        end
        
        function val = get.Velocity(this)
            val = this.angularVelocity;
        end
    end
    
    methods(Static, Access=public)
        function obj = Empty(simulationNumber)
            obj = AngularVelocityInBodyFrame(simulationNumber);
        end
    end
end