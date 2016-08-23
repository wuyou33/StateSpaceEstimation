function [ x, y ] = odeEuler(odefun, xspan, seed, step)
    % ODEEULER Euler ODE solver.
    %   [ x, y ] = odeEuler(odefun, xspan, seed, step).
    %    find details here: https://en.wikipedia.org/wiki/Euler_method
    %    uses EULER'S method to integrate an ODE.
    %
    %   INPUT:
    %       odefun  - pointer to ode function
    %       tspan   - [ti,tf] where ti and tf = initial and final values of independent variables
    %       y0      - initial value of dependent variable
    %       h       - step size
    %       p1, p2  - additional parameter used by dydt
    %
    %   OUTPUT:
    %       t - vector of independent variable
    %       y - vector of solution for dependent variable
    
    if nargin<4; error('at least 4 input arguments required'); end
    
    ti = xspan(1); 
    tf = xspan(2);
    
    if ~ (tf > ti); error('upper limit must be greater than lower limit'); end
    
    x = (ti:step:tf)';
    n = length(x);
    
    % if necessary, add an additional value of t so that range goes from t=ti to tf
    if x(n) < tf
        x(n+1) = tf;
        n = n+1;
        x(n)=tf;
    end
    
    y = cvecrep(seed, n);
    
    for i = 1:n-1
        y(:, i+1) = y(:, i) + odefun(x(i), y(:, i)) *( x(i+1)-x(i) );
    end
    
end
