classdef AccelerometerParam < handle
    %%
    properties( Access = private)
        simulationNumber;
        levelArm;
        angularAccelerationinBodyFrame;
        accelerometerScale;
        biasMu;
        biasSigma;
        noiseVar;
        bias;
    end 
    %%
    properties (Dependent)        
        LevelArm;
        AngularAccelerationinBodyFrame;
        AccelerometerScale;
        Bias;
        NoiseVar;
    end
    %%
    methods
        %%
        function obj = AccelerometerParam(simulationNumber, ...
                                          levelArm, ...
                                          angularAccelerationinBodyFrame, ...
                                          accelerometerScale, ...
                                          biasMu, ...
                                          biasSigma, ... 
                                          noiseVar)
                                      
            obj.simulationNumber = simulationNumber;
            obj.levelArm = levelArm;
            obj.angularAccelerationinBodyFrame = angularAccelerationinBodyFrame;
            obj.accelerometerScale = accelerometerScale;
            obj.biasMu = biasMu;
            obj.biasSigma = biasSigma;
            obj.noiseVar = noiseVar;       
        end
        %%
        function val = get.LevelArm(obj)
            val = obj.levelArm;
        end;
        %%
        function val = get.AngularAccelerationinBodyFrame(obj)
            val = obj.angularAccelerationinBodyFrame;
        end
        %% 
        function val = get.AccelerometerScale(obj)
            val = obj.accelerometerScale;
        end
        %%
        function val = get.Bias(obj)
            if (isempty(obj.bias))
                bmModel = bm(obj.biasMu, obj.biasSigma);
                bmModel.StartState = 0;                      
                obj.bias = simulate(bmModel, obj.simulationNumber);  
            end
            
            val = obj.bias;
        end
        %% 
        function val = get.NoiseVar(obj)
            val = obj.noiseVar;
        end
    end
end