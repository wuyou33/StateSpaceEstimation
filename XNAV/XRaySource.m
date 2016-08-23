classdef XRaySource < handle
    % Describe source for x-Ray navigation.
    % Source can be Pulsar or Quasar (with high stable EM emitting)
    properties (Access = public)
        name;
        period;             % sec
        intensity;          % ph·sec–1·cm-2
        raError;            % acrsec
        decError;           % acrsec
        gSource;
        gBackgr;
        galacticLon;        % [deg]
        galacticLat;        % [deg]
        distance;           % [km]
        normal;             % normal to x-ray source in solar system
        twoPiOnPeriod;      % 2*pi/period [1/sec]
    end
    
    properties (Dependent)
        Normal;
        TwoPiOnPeriod;
    end
    
    methods
        % constructor
        function obj = XRaySource(name, period, intensity, raError, decError, gSource, gBackgr, galacticLon, galacticLat, distance)
            % distance in kpc
            narginchk(10, 10);
            
            obj.name           = name;
            obj.period         = period;
            obj.intensity      = intensity;
            obj.raError        = raError;
            obj.decError       = decError;
            obj.gSource        = gSource;
            obj.gBackgr        = gBackgr;
            obj.galacticLon    = galacticLon;
            obj.galacticLat    = galacticLat;
            obj.distance       = distance*3.086e+16; % 1kpc is 3.086e+16 km
            obj.normal         = normFromLatLon(galacticLat, galacticLon);
            obj.twoPiOnPeriod  = 2*pi/period;
        end
        
        function val = get.Normal(this)
            val = this.normal;
        end
        
        function val = get.TwoPiOnPeriod(this)
            val = this.twoPiOnPeriod;
        end
    end
end
