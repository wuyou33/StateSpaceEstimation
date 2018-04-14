classdef InertialMeasurementUnit < handle
    % InertialMeasurementUnit. Describe inertial measurement unit (accelerometer and gyro).
    % Provides a measurement of acceleration (m / s^2) from three body axes accelerometer in body fixed frame
    % and measurement of angular velocity (redian per second) from three body axes gyro in body fixed frame
    
    properties( Access = private)
        accelerometerParams;
        accelerationInBodyFrame;
        gyroParams;
        angularVelocityInBodyFrame;
        timeData;
        acceleration;
        angularVelocity;
    end
    
    properties (Constant)
        g = 9.780327e-3; % [km / sec^2]
    end
    
    properties (Dependent)
        AngularVelocity;
        Acceleration;
        AccelerometerBias;
        GyroBias;
    end
    
    methods
        
        function obj = InertialMeasurementUnit(accelerometerParams, ...
                gyroParams, ...
                accelerationInBodyFrame, ...
                angularVelocityInBodyFrame, ...
                timeData)
            if ~isa(accelerometerParams, 'AccelerometerParam'); error('accelerometerParams should be instance of the AccelerometerParam'); end
            
            if ~isa(gyroParams, 'GyroParam'); error('gyroParams should be instance of the GyroParam'); end
            
            if ~isa(accelerationInBodyFrame, 'AccelerationInBodyFrame'); error('accelerationInBodyFrame should be intance of the AccelerationInBodyFrame'); end
            
            if ~isa(angularVelocityInBodyFrame, 'AngularVelocityInBodyFrame'); error('angularVelocityInBodyFrame should be intance of the AngularVelocityInBodyFrame'); end
            
            if ~isa(timeData, 'TimeExt'); error('timeData should be intance of the TimeExt'); end
            
            obj.accelerometerParams         = accelerometerParams;
            obj.accelerationInBodyFrame     = accelerationInBodyFrame;
            obj.gyroParams                  = gyroParams;
            obj.angularVelocityInBodyFrame  = angularVelocityInBodyFrame;
            obj.timeData                    = timeData;
        end
        
        function val = get.AngularVelocity(this)
            if isempty(this.angularVelocity);
                this.initializeAngularVelocity();
            end
            
            val = this.angularVelocity;
        end
        
        function val = get.Acceleration(this)
            if isempty(this.acceleration)
                this.initializeAcceleration();
            end
            
            val = this.acceleration;
        end
        
        function val = getAngularVelocity(this, index)
            if isempty(this.angularVelocity);
                this.initializeAngularVelocity();
            end
            
            val = this.angularVelocity(:, index);
        end
        
        function val = getAcceleartion(this, index)
            if isempty(this.acceleration)
                this.initializeAcceleration();
            end
            
            val = this.acceleration(:, index);
        end
        
        function val = get.AccelerometerBias(this)
            val = this.accelerometerParams.Bias;
        end
        
        function val = get.GyroBias(this)
            val = this.gyroParams.GyroBias;
        end
    end
    
    methods(Access = private)
        function initializeAcceleration(this)
            n = this.timeData.SimulationNumber;
            
            la_angAcc = column_vector_replicate(cross(this.accelerometerParams.LevelArm, this.accelerometerParams.AngularAccelerationinBodyFrame), n);
            
            w_w_la = cross(this.angularVelocityInBodyFrame.Velocity, ...
                cross(this.angularVelocityInBodyFrame.Velocity, ...
                column_vector_replicate(this.accelerometerParams.LevelArm, n)) ...
                );
            
            a  = this.accelerationInBodyFrame.Acceleration;
            sa = this.accelerometerParams.AccelerometerScale;
            ba = this.accelerometerParams.Bias;
            na = column_vector_replicate(this.accelerometerParams.NoiseVar, n) .* randn(3, n);
            
            this.acceleration = sa*(la_angAcc + w_w_la + a) + ba + na;
        end
        
        function initializeAngularVelocity(this)
            n  = this.timeData.SimulationNumber;
            w  = this.angularVelocityInBodyFrame.Velocity;
            a  = this.accelerationInBodyFrame.Acceleration;
            ga = this.gyroParams.GyroGSensitiveBias * a / this.g;
            sw = this.gyroParams.GyroScaleFactor;
            bw = this.gyroParams.GyroBias;
            nw = column_vector_replicate(this.gyroParams.GyroNoiseVar, n) .* randn(3, n);
            
            this.angularVelocity = ga + sw*w +  bw + nw;
        end
    end
    
end