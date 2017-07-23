classdef AstroBody
    % The class contains information about coordinates of solar system bodies like planets, sun, moons comets etc.
    % The class contain following
    % x     - position x coordinate [km];
    % y     - position y coordinate [km];
    % z     - position z coordinate [km];
    % vx    - velocity x coordinate [km/sec];
    % vy    - velocity y coordinate [km/sec];
    % vz    - velocity z coordinate [km/sec];
    % mass  - mass [kg].
    
    properties(SetAccess = private)
        X
        Y
        Z
        VX
        VY
        VZ
        Mass
    end
    
    methods
        function obj = AstroBody(r, v, mass)
            narginchk(3, 3);
            
            obj.X       = r(1);
            obj.Y       = r(2);
            obj.Z       = r(3);
            obj.VX      = v(1);
            obj.VY      = v(2);
            obj.VZ      = v(3);
            obj.Mass    = mass;
        end
    end
end
