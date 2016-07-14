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
    accelerometerParam = AccelerometerParam(initArgs.timeData, ...
        initArgs.levelArm, ... 
        initArgs.angularAccelerBodyFrame, ... 
        ones(3) + initArgs.accScale, ... 
        initArgs.accBiasMu, ... 
        initArgs.accBiasSigma, ... 
        initArgs.accNoiseVar...
    );

    gyroParam = GyroParam(initArgs.timeData, ...
        initArgs.gyroGSensitiveBias , ...
        ones(3) + initArgs.gyroScale, ...
        initArgs.gyroBiasMu, ...
        initArgs.gyroBiasSigma, ...
        initArgs.gyroNoiseVar ...        
    );

    inertialMeasurementUnit = InertialMeasurementUnit(accelerometerParam, ...
        gyroParam, ...
        initArgs.accelerationInBodyFrame, ...
        initArgs.angularVelocityInBodyFrame, ...
        initArgs.timeData ...
    );
       
    ins = InertialNavigationSystem(inertialMeasurementUnit, initArgs.timeData);
    
    if (initArgs.visualize)
        figure(); 
        
        subplot(3, 2, 1);
        plot2(initArgs.timeData.RelTime, inertialMeasurementUnit.AngularVelocity, 'angular velocity', {'x axis', 'y axis', 'z axis'}, 'angular velocity, rad/sec');

        subplot(3, 2, 3);
        plot2(initArgs.timeData.RelTime, initArgs.angularVelocityInBodyFrame.Velocity, 'true angular velocity', {'x axis', 'y axis', 'z axis'}, 'angular velocity, rad/sec');
            
        subplot(3, 2, 5);
        plot2(initArgs.timeData.RelTime, inertialMeasurementUnit.GyroBias, 'gyro bias', {'x axis', 'y axis', 'z axis'}, 'bias, rad/sec');
        
        subplot(3, 2, 2);
        plot2(initArgs.timeData.RelTime, 1e3*inertialMeasurementUnit.Acceleration, 'acceleration', {'x axis', 'y axis', 'z axis'}, 'acceleration, m/sec^2');
            
        subplot(3, 2, 4);
        plot2(initArgs.timeData.RelTime, 1e3*initArgs.accelerationInBodyFrame.Acceleration, 'true acceleration', {'x axis', 'y axis', 'z axis'}, 'acceleration, m/sec^2');            
        
        subplot(3, 2, 6);
        plot2(initArgs.timeData.RelTime, inertialMeasurementUnit.AccelerometerBias, 'acceleration bias', {'x axis', 'y axis', 'z axis'}, 'bias, km/sec^2');            
    end
end