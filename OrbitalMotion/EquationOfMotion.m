function dy = EquationOfMotion(time, state, acceleration, angularVelocity, tEpoch, sampleTime, startTime, gravityModel, mass)
    % EquationOfMotion. Solve motion equation for spaceship.
    %
    %   INPUT
    %        time                   time in [sec];
    %        state                  state space vector of spaceship (distance in [km], velocity in [km/s], quaternion);
    %        acceleration           acceleration of spaceship in [km/s^2];
    %        angularVelocity        angular velocity of spaceship in [rad/sec];
    %        tEpoch                 time Epoch;
    %        sampleTime             sample time required only if equation used in batch mode (batch mode used in imitator only) [sec];
    %        startTime              the start time of the solving of the whole problem as a whole. 
    %                                   required only if equation used in batch mode (batch mode used in imitator only) [sec]';
    %        gravityModel           Solar system gravity model (allow to calculate gravity accelaration in solar system;
    %        mass                   mass of the spacecraft [kg].
    %
    %   OUTPUT
    %        dy   increment (diff) of state space vector of spaceship (distance in [km], velocity in [km/s], quaternion).
    %
    EarthRadius = 6378.136; % [km] - Earth's equatorial radius
    MuE = 398600.4418;      % [km^3/s^2] - Earth gravity const
    J2 = .00108262575;
    J3 = -.000002533;
    J4 = -0.000001616;
    
    sample = round((time - startTime) / sampleTime) + 1;
    
    if isvector(acceleration)
        f = acceleration;
    else
        f = acceleration(:, sample);
    end
    
    if isvector(angularVelocity)
        w = angularVelocity;
    else
        w = angularVelocity(:, sample);
    end
    
    r = sqrt( state(1)^2 + state(2)^2 + state(3)^2 ); % distance from Earth center to spaceship center mass
    po = EarthRadius / r;
    
    sunInfluence = SunInfluence(tEpoch, state(1:3));
    moonInfluence = MoonInfluence(tEpoch, state(1:3));
    
    quaternion = state(7:10);
    rotatedAccelearation = quaternionRotation(quaternion, f, 2);
    
    gravAcc = gravityModel.EvalGravityAcceleration(sample, state(1:3), mass);
    
    dy(1:3, 1) = state(4:6);
    dy(4, 1) = gravAcc(1) -(MuE*state(1)/r^3)*(...
        J2*3/2*po^2*(1-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*(state(3)/r)^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4) ...
        ) ...
        + rotatedAccelearation(1) ...
        + sunInfluence(1) ...
        + moonInfluence(1);
    dy(5, 1) = gravAcc(2) -(MuE*state(2)/r^3)*(...
        J2*3/2*po^2*(1-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*(state(3)/r)^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4) ...
        ) ...
        + rotatedAccelearation(2) ...
        + sunInfluence(2) ...
        + moonInfluence(2);
    dy(6, 1) = gravAcc(3) -(MuE*state(3)/r^3)*(...
        J2*3/2*po^2*(3-5*(state(3)/r)^2) ...
        + J3*po^3*5/2*(3-7*po^2)*state(3)/r ...
        - J4*po^4*5/8*(3-42*(state(3)/r)^2+63*(state(3)/r)^4)...
        ) ...
        + rotatedAccelearation(3) ...
        + sunInfluence(3) ...
        + moonInfluence(3);
    
    dy(7:10, 1) = -0.5 * quaternionMultiply(quaternion, [0; w]);
end
