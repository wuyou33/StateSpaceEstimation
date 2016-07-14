classdef InertialNavigationSystem < handle
    
    properties (Access = private)
        imu;
        timeData;
    end
    
    methods
        function obj = InertialNavigationSystem(imu, timeData)
            if ~isa(imu, 'InertialMeasurementUnit')
                error('imu should be instance of the InertialMeasurementUnit');
            end
            
            obj.imu         = imu;
            obj.timeData    = timeData;
        end
        
        function state = simulate(this, initial, sample, tEpoch)
            if nargin == 2
                state = this.resolveBatch(initial);
            elseif nargin == 4
                state = this.resolve(initial, sample, tEpoch);
            else
                error('incorrect number of input argument. should be 2 or 4');
            end
        end
        
        function acceleration = getAcceleration(this, sample)
            acceleration = this.imu.getAcceleartion(sample);
        end;
        
        function angularVelocity = getAngularVelocity(this, sample)
            angularVelocity = this.imu.getAngularVelocity(sample);
        end;
    end
    
    methods (Access = private)
        function stateMatrix = resolveBatch(this, initial)
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = zeros(10, this.timeData.SimulationNumber);
            insState = initial;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            for i = 1 : num
                startBlock = (i-1)*blockSize + 1;
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = currentEpoch(this.timeData.JD, tMoonSun);                
                
                len = endBlock - startBlock + 1;
                if i == 1
                    stateMatrix(:, 1) = insState;
                end
                
                for j = startSample:len
                    insState = this.resolve(insState, j + startBlock - 1, tEpoch);
                    stateMatrix(:, j + startBlock - 1) = insState;
                end
                                
                tMoonSun = tMoonSun + this.timeData.RefreshSunMoonInfluenceTime;
                startSample = 1;
            end
        end
        
        function state = resolve(this, initial, sample, tEpoch)
            a  = this.imu.getAcceleartion(sample);
            w  = this.imu.getAngularVelocity(sample);
            dt = this.timeData.SampleTime;
            
            tEnd = this.timeData.Time(sample);
            timeSpan = [tEnd-dt, tEnd];
            
            odeFun = @(t,y) EquationOfMotion(t, y, a, w, tEpoch);
            [~, tmp] = ode45(odeFun, timeSpan, initial, odeset('MaxStep', dt));
            state = tmp(end, :)';
            state(7:10) = quaternionNormalize(state(7:10));
        end
    end
end
