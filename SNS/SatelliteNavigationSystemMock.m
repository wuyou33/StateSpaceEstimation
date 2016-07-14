classdef SatelliteNavigationSystemMock < SatelliteNavigationSystem
    properties (Access = private)
        simulatedState;
    end
    
    methods (Access = public)
        function this = SatelliteNavigationSystemMock(realState)
            simulationNumber = length(realState);
            this.simulatedState = realState + [1e-2*randn(simulationNumber, 3) 1e-4*randn(simulationNumber, 3) ];
        end
        
        function val = Simulate(this, sampleNumber)
            val = this.simulatedState(sampleNumber, :)';
        end
    end    
end

