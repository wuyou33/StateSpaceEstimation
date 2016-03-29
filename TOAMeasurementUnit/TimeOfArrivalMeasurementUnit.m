classdef TimeOfArrivalMeasurementUnit
        
    properties
        pulsar;
        observationEquation;
    end
    
    methods
        function obj = TimeOfArrivalMeasurementUnit(pulsar)
            obj.pulsar = pulsar;
            obj.observationEquation = @(state, pulsar) ...
            dot(state(1:3, :), ...
                repmat(pulsar.normalToSSB, ...
                    [1, ...
                     tern(length(state(1:3, :)) <= 3, 1, @() length(state(1:3, :)))]...
                    )...
                ) / velocitySpeed;
        end 
        
        function val = TimeOfArriaval(obj, state)
            val = obj.observationEquation(state);
        end
    end
    
end

