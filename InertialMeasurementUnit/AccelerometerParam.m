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
        sampleTime;
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
                                          noiseVar, ... 
                                          sampleTime)
                                      
            obj.simulationNumber = simulationNumber;
            obj.levelArm = levelArm;
            obj.angularAccelerationinBodyFrame = angularAccelerationinBodyFrame;
            obj.accelerometerScale = accelerometerScale;
            obj.biasMu = biasMu;
            obj.biasSigma = biasSigma;
            obj.noiseVar = noiseVar;  
            obj.sampleTime = sampleTime;
        end
        %%
        function val = get.LevelArm(this)
            val = this.levelArm;
        end;
        %%
        function val = get.AngularAccelerationinBodyFrame(this)
            val = this.angularAccelerationinBodyFrame;
        end
        %% 
        function val = get.AccelerometerScale(this)
            val = this.accelerometerScale;
        end
        %%
        function val = get.Bias(this)
            if (isempty(this.bias))
                % bmModel = bm(this.biasMu, this.biasSigma);
                % bmModel.StartState = 0;
                % this.bias = simulate(bmModel, this.simulationNumber);             
                dw = WienerProcess(this.biasMu, this.biasSigma);
                this.bias = dw.simulate(this.sampleTime, this.simulationNumber);
            end
            
            val = this.bias;
        end
        %% 
        function val = get.NoiseVar(this)
            val = this.noiseVar;
        end
    end
end