classdef TrajectoryPhaseSpaceSatelliteSimulator < handle
    properties (Access = private)
        initialState;
        accelerationInBodyFrame;
        angularVelocityInBodyFrame;
        timeTillCurrentEpoch;
    end
    
    methods 
        function obj = TrajectoryPhaseSpaceSatelliteSimulator(initialState, ...
                accelerationInBodyFrame, ...
                angularVelocityInBodyFrame, ...
                timeTillCurrentEpoch)
            
            if (~isa(accelerationInBodyFrame, 'AccelerationInBodyFrame'))
                error('accelerationInBodyFrame should be accelerometerParams of AccelerationInBodyFrame');
            end
            
            if (~isa(angularVelocityInBodyFrame, 'AngularVelocityInBodyFrame'))
                error('angularVelocityInBodyFrame should be accelerometerParams of AngularVelocityInBodyFrame');
            end
            
            obj.initialState = initialState;
            obj.accelerationInBodyFrame = accelerationInBodyFrame;
            obj.angularVelocityInBodyFrame = angularVelocityInBodyFrame;
            obj.timeTillCurrentEpoch = timeTillCurrentEpoch;
        end
        
        function res = Simulate(obj, time, sampleTime)
            [~, tmp] = ode113( @(t,y) EquationOfMotion(t, ...
                    y, ...
                    obj.accelerationInBodyFrame.Acceleration, ...
                    obj.angularVelocityInBodyFrame.Velocity, ...
                    obj.timeTillCurrentEpoch, ...
                    sampleTime ), ...
                time, ...
                obj.initialState ...
            );
            res = SatellitePhaseSpace(tmp', length(tmp));
        end
    end
    
end

