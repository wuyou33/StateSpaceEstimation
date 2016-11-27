function [ dy ] = equationOfMotionFreeFly( ~, state, tEpoch )
    % equationOfMotionFreeFly. Solve state dynamic equation for spaceship.
    %
    %   [ dy ] = equationOfMotionFreeFly( ~, state, tEpoch )
    %
    %   INPUT
    %        time       time in [sec];
    %        state      state space vector of spaceship (distance in [km], velocity in [km/s], quaternion);
    %        tEpoch     time Epoch.
    %
    %   OUTPUT
    %        dy     increment (diff) of state space vector of spaceship (distance in [km], velocity in [km/s]).
    %
    EarthRadius = 6378.136; % [km] - Earth's equatorial radius
    MuE = 398600.4418;      % [km^3/s^2] - Earth gravity const
    J2  = 0.00108262575;
    J3  = -0.000002533;
    J4  = -0.000001616;
    
    r = sqrt( state(1)^2 + state(2)^2 + state(3)^2 ); % distance from Earth center to spaceship center mass
    po = EarthRadius / r;
    
    dy = zeros(length(state), 1);
    
    sunInfluence = SunInfluence(tEpoch, state(1:3));
    moonInfluence = MoonInfluence(tEpoch, state(1:3));
    
    dy(1:3) = state(4:6);
    dy(4) = -(MuE*state(1)/r^3)*(...
        1 ...
        + J2*3/2*po^2*(1-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*(state(3)/r)^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4) ...
        ) ...
        + sunInfluence(1) ...
        + moonInfluence(1);
    dy(5) = -(MuE*state(2)/r^3)*(...
        1 ...
        + J2*3/2*po^2*(1-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*(state(3)/r)^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4) ...
        ) ...
        + sunInfluence(2) ...
        + moonInfluence(2);
    dy(6) = -(MuE*state(3)/r^3)*(...
        1 ...
        + J2*3/2*po^2*(3-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*po^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4)...
        ) ...
        + sunInfluence(3) ...
        + moonInfluence(3);
end
