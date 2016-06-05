function [varargout] = initInertialNavigationSystem(func, varargin)

    switch func  

        case 'init'
            model = init(varargin{1});
            varargout{1} = model;        

        otherwise        
            error(['Function ''' func ''' not supported.']);

    end
end

function ins = init(initArgs)
    accelerometerParam = AccelerometerParam(initArgs.simulationNumber, ... simulationNumber
        [0 0 0], ... levelArm
        [0 0 0], ... angularAccelerationinBodyFrame
        1*eye(3) + initArgs.accScale, ... accelerometerScale
        initArgs.accBiasMu, ... biasMu
        initArgs.accBiasSigma, ... biasSigma
        initArgs.accNoiseVar, ... noiseVar
        initArgs.sampleTime ... sampleTime
    );

    gyroParam = GyroParam(initArgs.simulationNumber, ... simulationNumber
        0*eye(3), ... gyroGSensitiveBias
        1*eye(3) + initArgs.gyroScale, ... gyroScaleFactor
        initArgs.gyroBiasMu, ... gyroBiasMu
        initArgs.gyroBiasSigma, ... gyroBiasSigma
        initArgs.gyroNoiseVar, ... gyroNoiseVar
        initArgs.sampleTime ... sampleTime
    );

    inertialMeasurementUnit = InertialMeasurementUnit(accelerometerParam, ...
        gyroParam, ...
        initArgs.accelerationInBodyFrame, ...
        initArgs.angularVelocityInBodyFrame, ...
        initArgs.simulationNumber ...
    );
       
    ins = InertialNavigationSystem(inertialMeasurementUnit, ...
        @(t, y, acceleration, angularVelocity, sampleTime) EquationOfMotion(t, y, acceleration, angularVelocity, initArgs.T_till_current_epoch, sampleTime ));       
    
    if (initArgs.visualize)
        figure(); 
            plot(initArgs.timeMinutes', inertialMeasurementUnit.AngularVelocity); 
            title('angular velocity'); 
            ylabel('angular velocity, rad/sec');
            xlabel('time, min');
            legend('x axis', 'y axis', 'z axis');
            grid on;

        figure(); 
            plot(initArgs.timeMinutes', initArgs.angularVelocityInBodyFrame.Velocity); 
            title('true angular velocity'); 
            ylabel('angular velocity, rad/sec');
            xlabel('time, min');
            legend('x axis', 'y axis', 'z axis');
            grid on;

        figure(); 
            plot(initArgs.timeMinutes', inertialMeasurementUnit.Acceleration); 
            ylabel('acceleration, km/sec^2');
            xlabel('time, min');
            legend('x axis', 'y axis', 'z axis');
            title('acceleration'); 
            grid on;

        figure(); 
            plot(initArgs.timeMinutes', initArgs.accelerationInBodyFrame.Acceleration); 
            ylabel('acceleration, km/sec^2');
            xlabel('time, min');
            legend('x axis', 'y axis', 'z axis');
            title('true acceleration'); 
            grid on;
    end
end