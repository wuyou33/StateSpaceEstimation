classdef InertialNavigationSystem < handle
        
    properties (Access = private)
        dynamicEquation;
        inertialMeasurementUnit;
    end
    
    methods
        function obj = InertialNavigationSystem(inertialMeasurementUnit, dynamicEquation)
            if (~isa(inertialMeasurementUnit, 'InertialMeasurementUnit'))
                error('inertialMeasurementUnit should be instance of the InertialMeasurementUnit');
            end
            
            obj.dynamicEquation = dynamicEquation;
            obj.inertialMeasurementUnit = inertialMeasurementUnit;
        end;
        
        function state = Simulate(this, initial, sample, sampleTime, currentTime)
            acceleration = this.inertialMeasurementUnit.getAcceleartion(sample);
            angularVelocity = this.inertialMeasurementUnit.getAngularVelocity(sample);
            [~, tmp] = ode113( @(t,y) this.dynamicEquation(t, y, acceleration, angularVelocity, sampleTime ), [currentTime - sampleTime, currentTime], initial );
            state = tmp(end, :)';
%             state = rungeKuttaFourOrderWithFixedStep( @(t,y) this.dynamicEquation(t, y, acceleration, angularVelocity, sampleTime ), initial, currentTime, sampleTime );
        end;
        
        function acceleration = GetAcceleration(this, sample)
            acceleration = this.inertialMeasurementUnit.getAcceleartion(sample);
        end;
        
        function angularVelocity = GetAngularVelocity(this, sample)
            angularVelocity = this.inertialMeasurementUnit.getAngularVelocity(sample);
        end;
    end
    
end

