classdef AccelerationInBodyFrame < handle
    %%
    properties (Access = private)
        acceleration;
        model;
        simulationNumber;
    end
    %%
    properties (Dependent)
        Acceleration;
    end
    
    methods
        %%
        function obj = AccelerationInBodyFrame(mu, sigma, simulationNumber)
            bmModel = bm(mu, sigma);
            bmModel.StartState = 0;
            obj.acceleration = simulate(bmModel, simulationNumber);     
            obj.model = bmModel;
            obj.simulationNumber = simulationNumber;
        end
        %%
        function val = get.Acceleration(obj)
            val = obj.acceleration;
        end
    end
end