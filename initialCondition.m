function [ condition ] = initialCondition( initialCov )
%%  Build initial condition for integrated INS and SNS navigation system

%%
    condition = sqrt(initialCov)*randn(22, 1); % sqrt(diag(initialCov))
    condition(7) = 1;
    condition(7:10) = quaternionNormalize(condition(7:10));
end

