classdef AngularVelocityInBodyFrame < handle
    %% Provide angular velocity caused, which represent angular rotation of body fixed frame
    % [rad / sec]
    %%
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
        function obj = AngularVelocityInBodyFrame(timeData, mu, sigma)
            %             bmModel = bm(mu, sigma);
            %             bmModel.StartState = 0;
            %
            %             this.angularVelocity = simulate(bmModel, simulationNumber);
            %             this.model = bmModel;
            obj.simulationNumber = timeData.SimulationNumber;
            
            if (nargin == 1)
                obj.dt    = NaN;
                obj.model = NaN;
                obj.angularVelocity = zeros(timeData.SimulationNumber, AngularVelocityInBodyFrame.Dimension);
            else
                obj.dt = timeData.SampleTime;                
                dw = WienerProcess(mu, sigma);
                obj.model = dw;
                obj.angularVelocity = dw.simulate(timeData.SampleTime, timeData.SimulationNumber);
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