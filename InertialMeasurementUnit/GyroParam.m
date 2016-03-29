classdef GyroParam < handle
    %%
    properties (Access = private)
        gyroGSensitiveBias;
        simulationNumber;        
        gyroScaleFactor;
        gyroBiasMu;
        gyroBiasSigma;
        gyroNoiseVar;
        bias;
    end
    %%
    properties (Dependent)
        GyroGSensitiveBias;
        GyroScaleFactor;
        GyroBias;
        GyroNoiseVar;
    end
    %%
    methods
        %%
        function this = GyroParam(simulationNumber, ...
                                 gyroGSensitiveBias, ...                                 
                                 gyroScaleFactor, ... 
                                 gyroBiasMu, ...
                                 gyroBiasSigma, ...
                                 gyroNoiseVar)            
            this.simulationNumber = simulationNumber;
            this.gyroGSensitiveBias = gyroGSensitiveBias;
            this.gyroScaleFactor = gyroScaleFactor;
            this.gyroBiasMu = gyroBiasMu;
            this.gyroBiasSigma = gyroBiasSigma;
            this.gyroNoiseVar = gyroNoiseVar;
        end
        %%
        function val = get.GyroGSensitiveBias(this)
            val = this.gyroGSensitiveBias;
        end
        %%
        function val = get.GyroScaleFactor(this)
            val = this.gyroScaleFactor;
        end
        %%
        function val = get.GyroBias(this)
            if (isempty(this.bias))
                bmModel = bm(this.gyroBiasMu, this.gyroBiasSigma);
                bmModel.StartState = 0;
                this.bias = simulate(bmModel, this.simulationNumber);
            end
            val = this.bias;
        end
        %%
        function val = get.GyroNoiseVar(this)
            val = this.gyroNoiseVar;
        end 
    end 
end 