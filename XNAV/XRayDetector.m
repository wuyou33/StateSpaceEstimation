classdef XRayDetector < handle
    % Describe photon detector for x-Ray navigation.
    % Source can be Pulsar or Quasar (with high stable EM emitting)
    
    properties(Access = private)
        xRaySources;
        detectorArea;
        timeBucket;
        backgroundPhotonRate;
        earthEphemeris;
        sunEphemeris;
        timeData;
        signal;
        spaceshipState;
        errorBudget;
    end
    
    properties (Access = public, Dependent)
        DetectorArea;
        TimeBucket;
        BackgroundPhotnRate;
    end
    
    methods(Access = public)
        function obj = XRayDetector(args)
            % Create a new instance of the XRayDetector
            %   xRaySources   - array of x-ray sources;
            %   detectorArea  - detector area (cm^2);
            %   timeBucket    - size of the bins used to count photons, ie total observed time (sec);
            %   backgroundPhotonRate  - is the background rate (including detector noise, the diffuse X-ray background,
            %               un-cancelled cosmic ray events and steady emission from the pulsar, bc) (photons/cm^2/sec);
            %   earthEphemeris - Earth ephemeris;
            %   sunEphemeris   - Sun ephemeris;
            %   spaceshipState - state (trajectory and velocity of spaceship during simulation);
            %
            narginchk(1, 1);
            
            obj.xRaySources            = args.xRaySources;
            obj.detectorArea           = args.detectorArea;
            obj.timeBucket             = args.timeBucket;
            obj.backgroundPhotonRate   = args.backgroundPhotonRate;
            obj.earthEphemeris         = args.earthEphemeris;
            obj.sunEphemeris           = args.sunEphemeris;
            obj.timeData               = args.timeData;
            obj.spaceshipState         = args.spaceshipState;
            obj.errorBudget            = args.errorBudget;
            obj.signal = [];
        end
        
        function observations = toa(this, visualize)
            if isempty(this.signal)
                this.evaluateTOA(visualize);
            end
            
            observations = this.signal;
        end
        
        function observation = getTOA(this, sample)
            if isempty(this.signal)
                this.evaluateTOA(0);
            end
            
            observation = this.signal(:, sample);
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
            res = this.backgroundPhotonRate;
        end
    end
    
    methods(Access = private)
        function evaluateTOA(this, visualize)
            % calculate difference time of arrival (DTOA)
            narginchk(2, 2);
            
            dtoa = calculateDiffToa(this.xRaySources, this.earthEphemeris, this.sunEphemeris, this.spaceshipState(1:3, :));
            this.signal = dtoa + this.generateNoise(size(dtoa, 2));
            
            if visualize
                this.visualize();
            end
        end
        
        function invPeriods = getInvPeriods(this)
            invPeriods = getInvPeriods(this.xRaySources);
        end
        
        function noise = generateNoise(this, capacity)
            if capacity <= 0; error('capacity should be positive integer'); end
            if capacity - fix(capacity) ~= 0; error('capacity should be positive integer'); end
            
            covariance  = xRayToaCovariance(this.xRaySources, this.detectorArea, this.timeBucket, this.backgroundPhotonRate, this.errorBudget);
            sourceCount = length(this.xRaySources);
            noise = chol(covariance, 'lower')*randn(sourceCount, capacity);
        end
        
        function [] = visualize(this)
            legends = {length(this.xRaySources)};
            figure();
            
            for i = 1:length(this.xRaySources)
                plot(this.timeData.Time, this.signal(i, :), 'LineWidth', 1);
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
