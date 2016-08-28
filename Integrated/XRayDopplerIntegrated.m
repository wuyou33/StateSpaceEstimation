classdef XRayDopplerIntegrated < handle
    %XRAYDOPPLERINTEGRATED Integrated X-Ray navigation system and Doppler navigation system.
    %   More information see in IET Radar Sonar Navig., 2011, Vol. 5, Iss. 9, pp. 1010-1017.
    
    properties (Access = private)
        timeData;
        xRayNavSystem;
        dopplerNavSystem;
        dimension = 6;
    end
    
    methods (Access = public)
        function obj = XRayDopplerIntegrated(timeData, xRayNavSystem, dopplerNavSystem)
            narginchk(3, 3);
            
            obj.timeData = timeData;
            obj.xRayNavSystem = xRayNavSystem;
            obj.dopplerNavSystem = dopplerNavSystem;
        end
        
        function stateMatrix = resolve(this, xRayEstimatorType, dopplerEstimatorType, xRayInitState, xRayInitCov, dopplerInitState, dopplerInitCov, report)
            narginchk(8, 8);
            
            stateMatrix = zeros(this.dimension, this.timeData.SimulationNumber);
            
            stateMatrix(:, 1) = dopplerInitState;
            
            dState = dopplerInitState;
            dCov = dopplerInitCov;            
            this.dopplerNavSystem.init(dopplerEstimatorType);
            dDecCov = this.dopplerNavSystem.updateFilterParams(dCov, dopplerEstimatorType);
            dSet = this.dopplerNavSystem.initParticleSet(dopplerEstimatorType, dState, dCov, dDecCov);
            
            rState = xRayInitState;
            rCov = xRayInitCov;
            this.xRayNavSystem.init(xRayEstimatorType);
            rDecCov = this.xRayNavSystem.updateFilterParams(rCov, xRayEstimatorType);
            rSet = this.xRayNavSystem.initParticleSet(xRayEstimatorType, rState, rCov, rDecCov);
            
            tMoonSun = this.timeData.StartSecond;
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                len = endBlock - startBlock + 1*(i == 1);
                
                for j = startSample:len
                    k = j + startBlock - 1;
                    
                    if report && mod((k / (this.timeData.SimulationNumber - 1))*100, 10) == 0
                        disp(['Completed: ', num2str((k / (this.timeData.SimulationNumber - 1)) * 100),' %' ]);
                    end
                    
                    [dState, dCov, dDecCov, dSet] = this.dopplerNavSystem.estimate(dState, dCov, dDecCov, dopplerEstimatorType, dSet, k, tEpoch);
                    [rState, rCov, rDecCov, rSet, ~] = this.xRayNavSystem.estimate(rState, rCov, rDecCov, xRayEstimatorType, rSet, k, tEpoch);
                                        
                    [state, ~, dState, dCov, rState, rCov] = federatedFilter(dState, dCov, rState, rCov);
                    stateMatrix(:, i) = state;
                    
                    stateMatrix(:, k) = state;
                end
                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
        end
    end
    
end
