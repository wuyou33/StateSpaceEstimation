classdef TrajectoryPhaseSpaceSatelliteSimulatorFree
    % TrajectoryPhaseSpaceSatelliteSimulatorFree. Simulate motion of spacecraft in simple spacecraft phase space (without attitude estimation).
    
    properties (Access = private)
        timeData;
        dimension = 6;
        mass;
    end
    
    methods
        function obj = TrajectoryPhaseSpaceSatelliteSimulatorFree(timeData, mass)
            narginchk(2, 2);
            
            if ~isa(timeData, 'TimeExt');
                error('timeData should be instance of the TimeExt');
            end
            
            obj.timeData = timeData;
            obj.mass = mass;
        end
        
        function signal = simulate(this, initialState, visualize)
            narginchk(2, 3);
            
            signal = zeros(this.dimension, this.timeData.SimulationNumber);
            
            m_fit_solar_system_gravity_model = memoize(@fit_solar_system_gravity_model);
            gravity_model = m_fit_solar_system_gravity_model(this.timeData.SampleTime, this.timeData.SimulationNumber);
            
            signal(1) = initialState;
            for i = 2:this.timeData.SimulationNumber
                c_time = this.timeData.Time(i);
                c_epoch = this.timeData.get_current_epoch(c_time);
                
                signal(i) = this.evaluate(signal(i-1), gravity_model, c_epoch, c_time);
            end
            
            if nargin == 3 && visualize
                signalWithRotation = [signal; ones(1, size(signal, 2)); zeros(3, size(signal, 2))];
                satellite_orbit_visualization(SatellitePhaseSpace(signalWithRotation, size(signalWithRotation, 2)));
            end
        end
        
        function signal = evaluate(this, initial, gravModel, tEpoch, time)
            narginchk(5, 5);
            
            odeFun = @(t, y) equation_of_motion_free_fly(t, y, tEpoch, gravModel, this.mass, this.timeData.SampleTime, this.timeData.StartSecond);
            [~, tmp] = ode_euler(odeFun, time, initial, this.timeData.SampleTime);
            
            if size(time, 2) == 2
                signal = tmp(end, :)';
            else
                signal = tmp';
            end
        end
    end
    
end
