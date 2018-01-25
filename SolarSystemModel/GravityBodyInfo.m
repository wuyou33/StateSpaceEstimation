classdef GravityBodyInfo
    % GravityBodyInfo
    % Contains information for astro objects (planet, moon, sun, comet or another massive object in solar system)
    % which allow to calculate gravitational forces between the object and any other object with the well-known mass and gravity.
    % <GravityBodyInfo> contains following parameters:
    %   mass            mass of the object [kg];
    %   position        position vector in ICRF-j2000 coordinate system [km].
    %
    properties (SetAccess = private)
        Mass
        Position
    end
    
    methods
        function obj = GravityBodyInfo(mass, position)
            narginchk(2, 2);
            
            if (size(position, 1) ~= 3)
                error('[ GravityBodyInfo : ctor ] position should be a array of column vector');
            end
            
            obj.Mass        = mass;
            obj.Position    = position;
        end
    end
end
