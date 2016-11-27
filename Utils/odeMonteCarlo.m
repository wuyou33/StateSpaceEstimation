function [ integ ] = odeMonteCarlo( x0, xn, odeFun )
    % odeMonteCarlo. Calculate Monte Carlo integration of an any function odeFunc between x0 and xn interval.
    %
    %   [ integ ] = odeMonteCarlo( xo, xn )
    %
    %   INPUT
    %       x0          start range of independent variable;
    %       xn          end of range of independent variable;
    %       odeFun      pointer to ode function.
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
            x = x0 + (xn - x0).*rand;
            ssum = ssum + odeFun(x);
        end
        N = i*10;
        integ = ssum / N;
    end
end
