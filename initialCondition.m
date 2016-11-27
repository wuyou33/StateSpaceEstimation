function [ condition ] = initialCondition( initialCov )
    %  initialCondition. Build initial condition for integrated INS and SNS navigation system.
    %
    %   [ condition ] = initialCondition( initialCov )
    %
    %   INPUT
    %       initialCov  covariation matrix of initial conditions.
    %
    %   OUTPUT
    %       condition   initial conditions.
    %
    condition = chol(initialCov, 'lower') * randn(22, 1); % sqrt(diag(initialCov))
    condition(7) = 1;
    condition(7:10) = quaternionNormalize(condition(7:10));
end
