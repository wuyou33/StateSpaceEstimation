classdef XRayDetector < handle
    % Describe photon detector for x-Ray navigation.
    % Source can be Pulsar or Quasar (with high stable EM emitting)
    
    properties(Access = private)
        xRaySources;
        detectorArea;
        timeBucket;
        backgroundPhotnRate;
        initialSpaceshipState;
        earthEphemeris;
        sunEphemeris;
        tEpoh
    end
    
    properties (Access = public, Dependent)
        DetectorArea;
        TimeBucket;
        BackgroundPhotnRate;
    end
    
    methods(Access = public)
        function obj = XRayDetector(xRaySources, detectorArea, timeBucket, backgroundPhotnRate, initialSpaceshipState, earthEphemeris, sunEphemeris, tEpoh)
            %% Create a new instance of the XRayDetector
            %   xRaySources   - array of x-ray sources;
            %   detectorArea  - detector area (cm^2);
            %   timeBucket    - size of the bins used to count photons, ie total observed time (sec);
            %   backgroundPhotnRate  - is the background rate (including detector noise, the diffuse X-ray background,
            %               un-cancelled cosmic ray events and steady emission from the pulsar, bc) (photons/cm^2/sec);
            %   initialSpaceshipState - initial state (orbit) of spaceship
            %   earthEphemeris - Earth ephemeris
            %   sunEphemeris   - Sun ephemeris
            %   tEpoh          - time of Epoh
            %%
            obj.xRaySources            = xRaySources;
            obj.detectorArea           = detectorArea;
            obj.timeBucket             = timeBucket;
            obj.backgroundPhotnRate    = backgroundPhotnRate;
            obj.initialSpaceshipState  = initialSpaceshipState;
            obj.earthEphemeris         = earthEphemeris;
            obj.sunEphemeris           = sunEphemeris;
            obj.tEpoh                  = tEpoh;
        end
        
        function signal = Simulate(this, time, visualize)
            if (nargin ~= 3); error('[ XRayDetector:Simulate ] incorrect number of input arguments'); end
            
            % calculate time of arrival (toa)
            simulator = TrajectoryPhaseSpaceSatelliteSimulatorFree(this.initialSpaceshipState, this.tEpoh);
            spaceshipState = simulator.Simulate(time);
            spaceshipTrajectory = spaceshipState(1:3, :)';
            diffToa = calculateDiffToa(this.xRaySources, this.earthEphemeris, this.sunEphemeris, spaceshipTrajectory);
            
            phase = diffToa2phase(this.getInvPeriods(), diffToa);
                        
            [capacity, ~] = size(phase);
            noise = this.generateNoise(capacity);
            
            signal = phase + noise;
            
            if (visualize)
                this.visualize(signal, time);
            end
        end
    end
    
    methods
        function res = get.DetectorArea(this)
            res = this.detectorArea;
        end
        
        function res = get.TimeBucket(this)
            res = this.timeBucket;
        end
        
        function res = get.BackgroundPhotnRate(this)
            res = this.backgroundPhotnRate;
        end
    end
    
    methods(Access = private)
        function invPeriods = getInvPeriods(this)
            dimension = length(this.xRaySources);            
            invPeriods   = zeros(1, dimension);
            for i = 1:dimension
                x = this.xRaySources(i);
                invPeriods(i) = x.TwoPiOnPeriod;
            end
        end
        
        function noise = generateNoise(this, capacity)
            if (capacity <= 0); error('capacity should be positive integer'); end
            if (capacity-fix(capacity) ~= 0); error('capacity should be positive integer'); end
            
            covariance  = xRayToaCovariance(this.xRaySources, this.detectorArea, this.timeBucket, this.backgroundPhotnRate);
            sourceCount = length(this.xRaySources);
            noise = randn(capacity, sourceCount)*diag(sqrt(covariance));
        end
        
        function [] = visualize(this, signal, time)     
            legends = {length(this.xRaySources)};
            figure();
            for i = 1:length(this.xRaySources)               
                plot(time, signal(:, i));
                hold on;
                legends{i} = this.xRaySources(i).name;
            end
            grid on;
            xlabel('time, sec');
            ylabel('phase, rad');
            legend(legends);
        end
    end
end