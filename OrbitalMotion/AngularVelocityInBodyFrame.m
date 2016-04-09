classdef AngularVelocityInBodyFrame < handle
    %%
    properties (Access = private)
        angularVelocity;
        model;
        simulationNumber;
        dt
    end
    %%
    properties (Dependent)
        Velocity;
    end
    %%
    methods
        function obj = AngularVelocityInBodyFrame(mu, sigma, simulationNumber, dt)
%             bmModel = bm(mu, sigma);
%             bmModel.StartState = 0;
%             
%             this.angularVelocity = simulate(bmModel, simulationNumber);
%             this.model = bmModel;
            obj.simulationNumber = simulationNumber;
            obj.dt = dt;
            
            dw = WienerProcess(mu, sigma);
            obj.model = dw;
            obj.angularVelocity = dw.simulate(dt, simulationNumber);
        end
        
        function val = get.Velocity(this)
            val = this.angularVelocity;
        end
    end
end