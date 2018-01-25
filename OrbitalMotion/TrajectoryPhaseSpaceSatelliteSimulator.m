classdef TrajectoryPhaseSpaceSatelliteSimulator < handle
    % TrajectoryPhaseSpaceSatelliteSimulator. Simulate motion of spacecraft in spacecraft phase space.
    
    properties (Access = private)
        initialState;
        accelerationInBodyFrame;
        angularVelocityInBodyFrame;
        timeData;
        mass;
        dimension = 10;
    end
    
    methods
        function obj = TrajectoryPhaseSpaceSatelliteSimulator(initialState, ...
                accelerationInBodyFrame, ...
                angularVelocityInBodyFrame, ...
                timeData, ...
                mass)
            
            narginchk(5, 5);
            
            if ~isa(accelerationInBodyFrame, 'AccelerationInBodyFrame')
                error('accelerationInBodyFrame should be instance of AccelerationInBodyFrame');
            end
            
            if ~isa(angularVelocityInBodyFrame, 'AngularVelocityInBodyFrame')
                error('angularVelocityInBodyFrame should be instance of AngularVelocityInBodyFrame');
            end
            
            if ~isa(timeData, 'TimeExt')
                error('timeData should be instance of the TimeExt');
            end
            
            obj.initialState = initialState;
            obj.accelerationInBodyFrame = accelerationInBodyFrame;
            obj.angularVelocityInBodyFrame = angularVelocityInBodyFrame;
            obj.timeData = timeData;
            obj.mass = mass;
        end
        
        function signal = simulate(this, visualize)
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = zeros(this.dimension, this.timeData.SimulationNumber);
            insState = this.initialState;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);
            gravModel = m_fitSolarSystemGravityModel(this.timeData.SampleTime, this.timeData.SimulationNumber);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                len = endBlock - startBlock + 1;
                if i == 1
                    stateMatrix(:, 1) = insState;
                end
                
                for j = startSample:len
                    insState = this.resolve(insState, j + startBlock - 1, tEpoch, gravModel);
                    stateMatrix(:, j + startBlock - 1) = insState;
                end
                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
            
            signal = SatellitePhaseSpace(stateMatrix, this.timeData.SimulationNumber);
            
            if nargin == 2 && visualize
                SatelliteOrbitVisualization(signal);
            end
        end
    end
    
    methods(Access = private)
        function state = resolve(this, initial, sample, tEpoch, gravityModel)
            a  = this.accelerationInBodyFrame.Acceleration(:, sample);
            w  = this.angularVelocityInBodyFrame.Velocity(:, sample);
            m  = this.mass;
            dt = this.timeData.SampleTime;
            t0 = this.timeData.StartSecond;
            
            tEnd = this.timeData.Time(sample);
            timeSpan = [tEnd-dt, tEnd];
            
            odeFun = @(t, y) EquationOfMotion(t, y, a, w, tEpoch, dt, t0, gravityModel, m);
            [~, tmp] = odeEuler(odeFun, timeSpan, initial, dt);
            
            state = tmp(end, :)';
            state(7:10) = quaternionNormalize(state(7:10));
        end
    end
    
end
