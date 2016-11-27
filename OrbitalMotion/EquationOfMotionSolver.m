function [ state ] = EquationOfMotionSolver(initialState, acceleration, angularVelocity, time, currentEpoch, sampleTime)
    % EquationOfMotionSolver. Solve DU which describe dynamic of spacecraft.
    %
    %   [ state ] = EquationOfMotionSolver(initialState, acceleration, angularVelocity,  time, T_till_current_epoch, sampleTime)
    %
    %   INPUT
    %       initialState        initial state;
    %       acceleration        all non-gravitational acceleration;
    %       angularVelocity     angular velocity;
    %       time                time span;
    %       currentEpoch        current Epoch;
    %       sampleTime          sample time.
    %
    %   OUTPUT
    %       state   calculated state vector of spacecraft.
    %
    [~, tmp] = ode113( @(t,y) EquationOfMotion(t, y, acceleration, angularVelocity, currentEpoch, sampleTime ), time, initialState );
    state = SatellitePhaseSpace(tmp');
end
