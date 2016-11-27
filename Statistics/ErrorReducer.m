classdef ErrorReducer < handle
    %% Provides methods to anlyze errors in estimated orbit parameters of spacecraft (estimated state space of satellity).
    
    properties (Access = private)
        iterations;
        reference;
        simulationNumber;
        iterationNumber;
        dimension;
    end
    
    methods
        function obj = ErrorReducer(iterations, reference, dimension, simulationNumber)
            obj.iterations = iterations;
            obj.reference = reference;
            obj.simulationNumber = simulationNumber;
            obj.iterationNumber = length(obj.iterations);
            obj.dimension = dimension;
        end
        
        function val = RMSD(obj)
            centredEstimate = zeros(obj.iterationNumber, obj.dimension);
            variance = zeros(obj.simulationNumber, obj.dimension);
            
            for i = 1:obj.simulationNumber
                for j = 2:obj.iterationNumber
                    centredEstimate(j,:) = obj.iterations(j).GetPhaseState(i)' - obj.reference.GetPhaseState(i)';
                end
                variance(i,:) = (1 / (obj.iterationNumber - 1)) * sum(centredEstimate.^2);
            end
            
            val = variance;
        end
    end
end
