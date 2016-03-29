classdef AngularVelocityInBodyFrame < handle
    %%
    properties (Access = private)
        angularVelocity;
        model;
        simulationNumber;
    end
    %%
    properties (Dependent)
        Velocity;
    end
    %%
    methods
        function this = AngularVelocityInBodyFrame(mu, sigma, simulationNumber)
            bmModel = bm(mu, sigma);
            bmModel.StartState = 0;
            
            this.angularVelocity = simulate(bmModel, simulationNumber);
            this.model = bmModel;
            this.simulationNumber = simulationNumber;
        end
        
        function val = get.Velocity(this)
            val = this.angularVelocity;
        end
    end
end