classdef InertialMeasurementUnit < handle
    % provides a measurement of acceleration (m / s^2) from three body axes accelerometer in body fixed frame
    % and measurement of angular velocity (redian per second) from three body axes gyro in body fixed frame
    
    %%
    properties( Access = private)
        accelerometerParams;
        accelerationInBodyFrame;
        gyroParams;
        angularVelocityInBodyFrame;
        simulationNumber;
        acceleration;
        angularVelocity;
    end
    
    %%
    properties (Constant)
        g = 9.780327;
    end
    
    %%
    properties (Dependent)
        AngularVelocity;
        Acceleration;
    end
    
    %%
    methods
        %%
        function obj = InertialMeasurementUnit(accelerometerParams, ...
                                               gyroParams, ...
                                               accelerationInBodyFrame, ...
                                               angularVelocityInBodyFrame, ...
                                               simulationNumber)
            if (~isa(accelerometerParams, 'AccelerometerParam'))
                error('params should be accelerometerParams of AccelerometerParam');
            end
            
            if (~isa(gyroParams, 'GyroParam'))
                error('params should be gyroParams of GyroParam');
            end
            
            if (~isa(accelerationInBodyFrame, 'AccelerationInBodyFrame'))
                error('params should be accelerationInBodyFrame of AccelerationInBodyFrame');
            end
            
            if (~isa(angularVelocityInBodyFrame, 'AngularVelocityInBodyFrame'))
                error('params should be angularVelocityInBodyFrame of AngularVelocityInBodyFrame');
            end
            
            obj.accelerometerParams = accelerometerParams;
            obj.accelerationInBodyFrame = accelerationInBodyFrame;
            obj.gyroParams = gyroParams;
            obj.angularVelocityInBodyFrame = angularVelocityInBodyFrame;
            obj.simulationNumber = simulationNumber;
            
            obj.acceleration = (...
                        repmat(...
                            cross(obj.accelerometerParams.LevelArm, obj.accelerometerParams.AngularAccelerationinBodyFrame), ...
                                obj.simulationNumber, ...
                                1) ...
                        + cross(obj.angularVelocityInBodyFrame.Velocity, ...
                            cross(obj.angularVelocityInBodyFrame.Velocity, ...
                                repmat(obj.accelerometerParams.AngularAccelerationinBodyFrame, obj.simulationNumber, 1 ) ...
                                ) ...
                            ) ... 
                        + obj.accelerationInBodyFrame.Acceleration ...
                    ) * obj.accelerometerParams.AccelerometerScale ... 
                    + obj.accelerometerParams.Bias ...
                    + (repmat(obj.accelerometerParams.NoiseVar, 1, obj.simulationNumber) .* randn(3, obj.simulationNumber))';
                
            obj.angularVelocity = obj.accelerationInBodyFrame.Acceleration / obj.g * obj.gyroParams.GyroGSensitiveBias ...
                + obj.angularVelocityInBodyFrame.Velocity * obj.gyroParams.GyroScaleFactor ...
                + obj.gyroParams.GyroBias ...
                + (repmat(obj.gyroParams.GyroNoiseVar, 1, obj.simulationNumber) .* randn(3, obj.simulationNumber))';
        end
        %%
        function val = get.AngularVelocity(this)
            val = this.angularVelocity;
        end
        %%
        function val = get.Acceleration(this)
            val = this.acceleration;
        end
        %%
        function val = getAngularVelocity(this, index)
            val = this.angularVelocity(index, :);
        end
        %%
        function val = getAcceleartion(this, index)
            val = this.acceleration(index, :);
        end
    end
end