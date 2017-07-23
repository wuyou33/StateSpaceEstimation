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
            
            if ~isa(timeData, 'TimeExt'); error('timeData should be instance of the TimeExt'); end
            
            obj.timeData = timeData;
            obj.mass = mass;
        end
        
        function signal = simulate(this, initialState, visualize)
            narginchk(2, 3);
            
            tMoonSun = this.timeData.StartSecond;
            signal = zeros(this.dimension, this.timeData.SimulationNumber);
            initial = initialState;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            
            m_fitSolarSystemGravityModel = memoize(@fitSolarSystemGravityModel);                        
            gravModel = m_fitSolarSystemGravityModel(this.timeData.SampleTime, this.timeData.SimulationNumber);
            
            for i = 1:num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                cTime = this.timeData.Time(startBlock:endBlock);
                cEpoch = currentEpoch(this.timeData.JD, tMoonSun);
                
                signal(:, startBlock:endBlock) = this.evaluate(initial, gravModel, cEpoch, cTime);
                
                initial = signal(:, endBlock);
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
            end
            
            if nargin == 3 && visualize
                signalWithRotation = [signal; ones(1, size(signal, 2)); zeros(3, size(signal, 2))];
                SatelliteOrbitVisualization(SatellitePhaseSpace(signalWithRotation, size(signalWithRotation, 2)));
            end
        end
        
        function signal = evaluate(this, initial, gravModel, tEpoch, time)
            narginchk(5, 5);
            
            odeFun = @(t, y) equationOfMotionFreeFly(t, y, tEpoch, gravModel, this.mass, this.timeData.SampleTime, this.timeData.StartSecond);
            [~, tmp] = ode45(odeFun, time, initial, odeset('MaxStep', this.timeData.SampleTime));
            
            if size(time, 2) == 2
                signal = tmp(end, :)';
            else
                signal = tmp';
            end
        end
    end
    
end
