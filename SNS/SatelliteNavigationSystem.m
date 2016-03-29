classdef (Abstract) SatelliteNavigationSystem < handle          
    methods (Abstract, Access = public)
        val = Simulate(obj, sampleNumber);
    end    
end

