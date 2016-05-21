classdef TrajectoryPhaseSpaceSatelliteSimulatorFree
    properties (Access = private)
        initialState;
        timeTillCurrentEpoch;
    end
    
    methods 
        function obj = TrajectoryPhaseSpaceSatelliteSimulatorFree(initialState, timeTillCurrentEpoch)            
            obj.initialState = initialState;            
            obj.timeTillCurrentEpoch = timeTillCurrentEpoch;
        end
        
        function signal = Simulate(this, time)
            [~, tmp] = ode45( @(t,y) equationOfMotionFreeFly(t, ...
                    y, ...
                    this.timeTillCurrentEpoch), ...
                time, ...
                this.initialState ...
            );
            signal = tmp';
        end
    end
    
end