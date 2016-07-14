classdef AccelerometerParam < handle
    
    properties( Access = private)
        timeData;
        levelArm;
        angularAccelerationinBodyFrame;
        accelerometerScale;
        biasMu;
        biasSigma;
        noiseVar;
        bias;
    end 
    
    properties (Dependent)        
        LevelArm;
        AngularAccelerationinBodyFrame;
        AccelerometerScale;
        Bias;
        NoiseVar;
    end
    
    methods
        
        function obj = AccelerometerParam(timeData, ...
                                          levelArm, ...
                                          angularAccelerationinBodyFrame, ...
                                          accelerometerScale, ...
                                          biasMu, ...
                                          biasSigma, ... 
                                          noiseVar)
                                      
            obj.timeData = timeData;
            obj.levelArm = levelArm;
            obj.angularAccelerationinBodyFrame = angularAccelerationinBodyFrame;
            obj.accelerometerScale = accelerometerScale;
            obj.biasMu = biasMu;
            obj.biasSigma = biasSigma;
            obj.noiseVar = noiseVar;
        end
        
        function val = get.LevelArm(this)
            val = this.levelArm;
        end;
        
        function val = get.AngularAccelerationinBodyFrame(this)
            val = this.angularAccelerationinBodyFrame;
        end
        
        function val = get.AccelerometerScale(this)
            val = this.accelerometerScale;
        end
        
        function val = get.Bias(this)
            if (isempty(this.bias))             
                dw = WienerProcess(this.biasMu, this.biasSigma);
                this.bias = dw.simulate(this.timeData.SampleTime, this.timeData.SimulationNumber);
            end
            
            val = this.bias;
        end
        
        function val = get.NoiseVar(this)
            val = this.noiseVar;
        end
    end
end