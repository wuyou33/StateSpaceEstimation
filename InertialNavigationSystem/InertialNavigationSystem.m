classdef InertialNavigationSystem < handle
    % InertialNavigationSystem. Inertial navigation system. Allow to estimate trajectory of spacecraft in full spacecraft phase space.
    
    properties (Access = private)
        imu;
        timeData;
        dimension = 10;
        mass;
        gravityModel;
    end
    
    properties (Dependent)
        SimulationNumber;
    end
    
    methods
        function obj = InertialNavigationSystem(imu, timeData, gravityModel, mass)
            narginchk(4, 4);
            
            if ~isa(imu, 'InertialMeasurementUnit')
                error('imu should be instance of the InertialMeasurementUnit');
            end
            
            obj.imu             = imu;
            obj.mass            = mass; % [kg]
            obj.timeData        = timeData;
            obj.gravityModel    = gravityModel;
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
        
        function stateMatrix = evaluate(this, initial, timeStart, timeEnd)
            narginchk(4, 4);
            
            startSample = this.timeData.evalSampleFromTime(timeStart);
            fromSample  = startSample + 1;
            endSample   = this.timeData.evalSampleFromTime(timeEnd);
            
            state = initial;
            if fromSample == endSample
                localTime = [0; timeStart];
            else
                localTime = timeStart : this.timeData.SampleTime : timeEnd - this.timeData.SampleTime;
            end
            
            stateMatrix = zeros(this.dimension, endSample - startSample);
            for i = fromSample : endSample
                timeAtEndOfSample = localTime(i - startSample);
                tMoonSunDelta = floor(timeAtEndOfSample / this.timeData.RefreshSunMoonInfluenceTime) * this.timeData.RefreshSunMoonInfluenceTime;
                tMoonSun = this.timeData.StartSecond + tMoonSunDelta;
                
                state = this.resolve(state, i, current_epoch(this.timeData.JD, tMoonSun));
                stateMatrix(:, i - startSample) = state;
            end
        end
        
        function acceleration = getAcceleration(this, sample)
            acceleration = this.imu.getAcceleartion(sample);
        end
        
        function angularVelocity = getAngularVelocity(this, sample)
            angularVelocity = this.imu.getAngularVelocity(sample);
        end
        
        function val = get.SimulationNumber(this)
            val = this.timeData.SimulationNumber;
        end
        
    end
    
    methods (Access = private)
        function stateMatrix = resolveBatch(this, initial)
            tMoonSun = this.timeData.StartSecond;
            stateMatrix = zeros(this.dimension, this.timeData.SimulationNumber);
            insState = initial;
            
            num = ceil(this.timeData.TotalSeconds / this.timeData.RefreshSunMoonInfluenceTime);
            blockSize = ceil(this.timeData.SimulationNumber / num);
            startSample = 2;
            
            for i = 1 : num
                startBlock = (i-1)*blockSize + 1*(i == 1);
                endBlock   = min(i*blockSize, this.timeData.SimulationNumber);
                
                tEpoch = current_epoch(this.timeData.JD, tMoonSun);
                
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
            m  = this.mass;
            dt = this.timeData.SampleTime;
            t0 = this.timeData.StartSecond;
            
            tEnd = this.timeData.Time(sample);
            timeSpan = [tEnd-dt, tEnd];
            
            odeFun = @(t,y) EquationOfMotion(t, y, a, w, tEpoch, dt, t0, this.gravityModel, m);
            [~, tmp] = ode_euler(odeFun, timeSpan, initial, dt);
            % [~, tmp] = ode45(odeFun, timeSpan, initial, odeset('MaxStep', dt));
            state = tmp(end, :)';
            state(7:10) = quaternion_normalize(state(7:10));
        end
    end
    
end
