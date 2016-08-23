classdef TrajectoryPhaseSpaceSatelliteSimulatorFree
    properties (Access = private)
        timeData;
        dimension = 6;
    end
    
    methods
        function obj = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeData)
            narginchk(1, 1);
            
            if ~isa(timeData, 'TimeExt'); error('timeData should be instance of the TimeExt'); end
            
            obj.timeData = timeData;
        end
        
        function signal = simulate(this, initialState)
            narginchk(2, 2);
            
            tMoonSun = this.timeData.StartSecond;
            signal = zeros(this.dimension, this.timeData.SimulationNumber);
            initial = initialState;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                signal(:, startBlock:endBlock) = this.evaluate(initial, currentEpoch(this.timeData.JD, tMoonSun), this.timeData.Time(startBlock:endBlock));
                
                initial = signal(:, endBlock);
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
            end
        end
        
        function signal = evaluate(this, initial, tEpoch, time)
            odeFun = @(t, y) equationOfMotionFreeFly(t, y, tEpoch);
            [~, tmp] = ode45(odeFun, time, initial, odeset('MaxStep', this.timeData.SampleTime));
            
            if size(time, 2) == 2
                signal = tmp(end, :)';
            else
                signal = tmp';
            end
        end
    end
    
end
