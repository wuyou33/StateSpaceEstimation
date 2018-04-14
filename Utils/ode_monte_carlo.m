function [ integ ] = ode_monte_carlo( x_init, xn, ode_fun )
    % ode_monte_carlo. Calculate Monte Carlo integration of an any function odeFunc between x_init and xn interval.
    %
    %   [ integ ] = ode_monte_carlo( xo, xn )
    %
    %   INPUT
    %       x_init      start range of independent variable;
    %       xn          end of range of independent variable;
    %       ode_fun     pointer to ode function.
    %
    %   OUTPUT
    %       integ   calculated result of integral.
    %
    % generating random number with the sequence reproducible
    rng('state', 0);
    
    ssum = 0;
    
    % calculate the function of a total of 1000 times
    % printing an estimate of the integral after every 10 iterations
    for i = 1:1000
        for j = 1:10
            x = x_init + (xn - x_init).*rand;
            ssum = ssum + ode_fun(x);
        end
        N = i*10;
        integ = ssum / N;
    end
end
