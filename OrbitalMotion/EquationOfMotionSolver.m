function [state] = EquationOfMotionSolver(initialState, acceleration, angularVelocity,  time, T_till_current_epoch, sampleTime)
    [~, tmp] = ode113( @(t,y) EquationOfMotion(t, y, acceleration, angularVelocity, T_till_current_epoch, sampleTime ), time, initialState );
    state = SatellitePhaseSpace(tmp');
end