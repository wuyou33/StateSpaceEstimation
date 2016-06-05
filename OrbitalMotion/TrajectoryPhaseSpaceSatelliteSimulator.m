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
        
        function signal = Simulate(this, time, sampleTime, visualize)
            [~, tmp] = ode113( @(t,y) EquationOfMotion(t, ...
                    y, ...
                    this.accelerationInBodyFrame.Acceleration, ...
                    this.angularVelocityInBodyFrame.Velocity, ...
                    this.timeTillCurrentEpoch, ...
                    sampleTime ), ...
                time, ...
                this.initialState ...
            );
            signal = SatellitePhaseSpace(tmp', length(tmp));
            
            if nargin == 4 && visualize == 1
                SatelliteOrbitVisualization(signal);
            end
        end
    end
    
end