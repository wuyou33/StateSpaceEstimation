function [ points, weights ] = evalFifthDegreeCubatureRule(n)
    % evalFifthDegreeCubatureRule. Evaluate cubature points and weights for 5-th degree Cubature Rule.
    %
    %   INPUT
    %       n - state space dimension;
    %
    %   OUPUT
    %       set     - evaluated cubature points;
    %       weights - corresponding cubature weights;
    %
    narginchk(1, 1);
    
    num = 2*n^2 + 1;
    points = zeros(n, num);
    weights = zeros(1, num);
    
    points(:, 1) = zeros(n, 1);
    weights(1) = 2 / (n + 2);
    
    for i = 1 : n
        points(i, i+1) = sqrt(n/2 + 1);
        weights(i+1) = (4 - n) / (2*(n + 2)^2);
        
        points(i, i + n + 1)  = -sqrt(n/2 + 1);
        weights(i + n + 1) = (4 - n) / (2*(n + 2)^2);
    end
    
    count = 2*n + 1;
    
    for i = 1 : n
        for j = (i + 1) : n
            count = count + 1;
            points(i, count) = sqrt(n/4 + 1/2);
            points(j, count) = sqrt(n/4 + 1/2);
            weights(count) = 1 / ( (n + 2)^2 );
            
            count = count + 1;
            points(i, count) = sqrt(n/4 + 1/2);
            points(j, count) = -sqrt(n/4 + 1/2);
            weights(count) = 1/((n+2)^2);
            
            count = count + 1;
            points(i, count) = -sqrt(n/4 + 1/2);
            points(j, count) = sqrt(n/4 + 1/2);
            weights(count) = 1 / ( (n + 2)^2 );
            
            count = count + 1;
            points(i, count) = -sqrt(n/4 + 1/2);
            points(j, count) = -sqrt(n/4 + 1/2);
            weights(count) = 1 / ( (n + 2)^2 );
        end
    end
    
    points = points*sqrt(2);
end
