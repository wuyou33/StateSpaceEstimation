classdef AccelerationInBodyFrame < handle
    %%
    properties (Access = private)
        acceleration;
        model;
        simulationNumber;
        dt
    end
    %%
    properties (Dependent)
        Acceleration;
    end
    
    methods
        %%
        function obj = AccelerationInBodyFrame(mu, sigma, simulationNumber, dt)
%             bmModel = bm(mu, sigma);
%             bmModel.StartState = 0;
%             obj.acceleration = simulate(bmModel, simulationNumber);     
%             obj.model = bmModel;
            obj.simulationNumber = simulationNumber;
            obj.dt = dt;
            
            dw = WienerProcess(mu, sigma);
            obj.model = dw;
            obj.acceleration = dw.simulate(dt, simulationNumber);
        end
        %%
        function val = get.Acceleration(obj)
            val = obj.acceleration;
        end
    end
end