classdef DopplerMeasurementUnit < handle
    %DOPPLERMEASUREMENTUNIT Doppler  measurement unit.
    %   The Doppler navigation utilises the Doppler shift to calculate the relative velocity between a fixed position and a spacecraft.
    %   This method usually adopts one of two Doppler shifts as the measurement:
    %   1) the Doppler shift caused by the relative motion from a spacecraft to a ground station  and
    %   2) the Doppler shift caused by the relative motion between a spacecraft and the Sun
    %   Implemented according to "IET Radar Sonar Navig., 2011, Vol. 5, Iss. 9, pp. 1010-1017"
    
    properties (Access = private)
        earthEphemeris;
        sunEphemeris;
        starshipState;
        timeData;
        variance;
        signal;
    end
    
    properties (Dependent)
        Variance;
    end
    
    methods (Access = public)
        function obj = DopplerMeasurementUnit(earthEphemeris, sunEphemeris, starshipState, timeData, variance)
            narginchk(5, 5);
            
            obj.earthEphemeris  = earthEphemeris;
            obj.sunEphemeris    = sunEphemeris;
            obj.starshipState   = starshipState;
            obj.timeData        = timeData;
            obj.variance        = variance;
        end
        
        function observation = simulate(this, visualize)
            narginchk(2, 2);
            
            if isempty(this.signal)
                this.resolve(visualize);
            end
            
            observation = this.signal;
        end
        
        function observation = dopplerShift(this, sample)
            narginchk(2, 2);
            
            if isempty(this.signal)
                this.resolve();
            end
            
            observation = this.signal(:, sample);
        end
    end
    
    methods
        function val = get.Variance(this)
            val = this.variance;
        end
    end
    
    methods (Access = private)
        function resolve(this, visualize)
            narginchk(1, 2);
            
            if nargin == 1
                visualize = 0;
            end
            
            this.signal = dopplerShift(this.starshipState, this.earthEphemeris, this.sunEphemeris) + ...
                sqrt(this.variance) * randn(1, this.timeData.SimulationNumber);
            
            if visualize
                this.visualize();
            end
        end
        
        function visualize(this)
            narginchk(1, 1);
            
            plot(this.timeData.Time, this.signal*1e3, 'LineWidth', 1);
            hold on;
            grid on;
            xlabel('time, sec');
            ylabel('radial speed, m / sec');
            legend('radial speed');
        end
    end
    
end
