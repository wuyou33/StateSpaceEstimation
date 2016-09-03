classdef GyroParam < handle
    
    properties (Access = private)
        gyroGSensitiveBias;
        gyroScaleFactor;
        gyroBiasMu;
        gyroBiasSigma;
        gyroNoiseVar;
        bias;
        timeData;
    end
    
    properties (Dependent)
        GyroGSensitiveBias;
        GyroScaleFactor;
        GyroBias;
        GyroNoiseVar;
    end
    
    methods
        
        function this = GyroParam(timeData, ...
                gyroGSensitiveBias, ...
                gyroScaleFactor, ...
                gyroBiasMu, ...
                gyroBiasSigma, ...
                gyroNoiseVar)
            this.timeData = timeData;
            this.gyroGSensitiveBias = gyroGSensitiveBias;
            this.gyroScaleFactor = gyroScaleFactor;
            this.gyroBiasMu = gyroBiasMu;
            this.gyroBiasSigma = gyroBiasSigma;
            this.gyroNoiseVar = gyroNoiseVar;
        end
        
        function val = get.GyroGSensitiveBias(this)
            val = this.gyroGSensitiveBias;
        end
        
        function val = get.GyroScaleFactor(this)
            val = this.gyroScaleFactor;
        end
        
        function val = get.GyroBias(this)
            if isempty(this.bias)
                dw = WienerProcess(this.gyroBiasMu, this.gyroBiasSigma);
                this.bias = dw.simulate(this.timeData.SampleTime, this.timeData.SimulationNumber);
            end
            
            val = this.bias;
        end
        
        function val = get.GyroNoiseVar(this)
            val = this.gyroNoiseVar;
        end
    end
end
