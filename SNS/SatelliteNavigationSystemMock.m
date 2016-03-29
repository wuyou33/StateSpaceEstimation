classdef SatelliteNavigationSystemMock < SatelliteNavigationSystem
    properties (Access = private)
        simulatedState;
    end
    
    methods (Access = public)
        function this = SatelliteNavigationSystemMock(realState)
            simulationNumber = length(realState);
            this.simulatedState = realState + [10*randn(simulationNumber, 3) 0.01*randn(simulationNumber, 3) ];
        end
        
        function val = Simulate(this, sampleNumber)
            val = this.simulatedState(sampleNumber, :)';
        end
    end    
end

