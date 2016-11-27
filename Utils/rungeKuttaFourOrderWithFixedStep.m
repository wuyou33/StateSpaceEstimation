function [ out ] = rungeKuttaFourOrderWithFixedStep( odeFun, y, x, span )
    % rungeKuttaFourOrderWithFixedStep. Runge-Kutta ODE solver (used fixed step size).
    %
    %   Here y_{n+1} is the RK4 approximation of y(t_{n+1}), and the next value (y_{n+1})
    %   is determined by the present value (y_n) plus the weighted average of four increments,
    %   where each increment is the product of the size of the interval, span,
    %   and an estimated slope specified by function f on the right-hand side of the differential equation.
    %
    %   k_1 is the increment based on the slope at the beginning of the interval, using dy(x, y)
    %   k_2 is the increment based on the slope at the midpoint of the interval, using dy(x+0.5*span, y+0.5*span*k_1)
    %   k_3 is again the increment based on the slope at the midpoint, but now using dy((x+0.5*span),(y+0.5*span*k_2));
    %   k_4 is the increment based on the slope at the end of the interval, using dy((x+span),(y+k_3*span));
    %
    %   [ out ] = rungeKuttaFourOrderWithFixedStep( odeFun, y, x, span )
    %
    %   INPUT
    %       odefun  	pointer to ode function;
    %       y           initial value of dependent variable;
    %       x           value of independent variable;
    %       span        [ti,tf] where ti and tf = initial and final values of independent variables.
    %
    %   OUTPUT
    %       out 	vector of solution for dependent variable.
    %
    k_1 = odeFun(x, y);
    k_2 = odeFun(x+0.5*span, y+0.5*span*k_1);
    k_3 = odeFun((x+0.5*span), (y+0.5*span*k_2));
    k_4 = odeFun((x+span), (y+k_3*span));
    
    out = y + (span/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
end
