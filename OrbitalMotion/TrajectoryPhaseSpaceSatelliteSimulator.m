classdef TrajectoryPhaseSpaceSatelliteSimulator < handle
    properties (Access = private)
        initialState;
        accelerationInBodyFrame;
        angularVelocityInBodyFrame;
        timeData;
        dimension = 10;
    end
    
    methods
        function obj = TrajectoryPhaseSpaceSatelliteSimulator(initialState, ...
                accelerationInBodyFrame, ...
                angularVelocityInBodyFrame, ...
                timeData)
            
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
        end
        
        function signal = simulate(this, visualize)
            tMoonSun = this.timeData.StartSecond;            
            stateVector = zeros(this.dimension, this.timeData.SimulationNumber);            
            initial = this.initialState;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                time = this.timeData.Time(startBlock:endBlock);
                
                a  = this.accelerationInBodyFrame.Acceleration;
                w  = this.angularVelocityInBodyFrame.Velocity;
                dt = this.timeData.SampleTime;
                t0 = this.timeData.StartSecond;
                
                odeFun = @(t, y) EquationOfMotion(t, y, a, w, tEpoch, dt, t0);                
                [~, tmp] = ode45(odeFun, time, initial, odeset('MaxStep', dt));
                                
                stateVector(:, startBlock:endBlock) = tmp';
                initial = stateVector(:, endBlock);
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
            end
            
            signal = SatellitePhaseSpace(stateVector, this.timeData.SimulationNumber);
            
            if nargin == 2 && visualize == 1
                SatelliteOrbitVisualization(signal);
            end
        end
    end
    
end