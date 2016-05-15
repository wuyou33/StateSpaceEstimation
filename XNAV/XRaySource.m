classdef XRaySource
    % Describe source for x-Ray navigation.
    % Source can be Pulsar or Quasar (with high stable EM emitting)        
    properties (Access = public)
        name;
        period;             %msec
        intensity;          % ph·sec–1·cm-2
        raError;            % acrsec
        decError;           % acrsec
        gSource;
        gBackgr;
        galacticLon;        % deg
        galacticLat;        % deg
        distance;           % km
    end
    
    properties (Dependent)
        Normal;
    end
    
    methods
        %% constructor
        function obj = XRaySource(name, period, intensity, raError, decError, gSource, gBackgr, galacticLon, galacticLat, distance)
            % distance in kpc
            if (nargin ~= 10)
                error('[XRaySource] not enough input args');
            end
            
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
        end
        
        function val = get.Normal(this)
            val = normFromLatLon(this.galacticLat, this.galacticLon);
        end
    end    
end