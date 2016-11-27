classdef (Abstract) SatelliteNavigationSystem < handle
    % Define interface for Satellite Navigation System, which can be used in integrated navigation system.
    
    methods (Abstract, Access = public)
        val = Simulate(obj, sampleNumber);
    end
    
end
