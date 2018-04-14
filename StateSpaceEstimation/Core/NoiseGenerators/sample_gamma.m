function [ noise ] = sample_gamma( noiseDataSet, count )
    % sample_gamma. Generate N (count) samples of a noise source specified by the noiseDataSet data structure (gamma stochastic processes).
    %
    %   [ noise ] = sample_gamma( noiseDataSet, count )
    %
    %   INPUT
    %       noiseDataSet    structure which fully describe stochastic process;
    %       count           count of requested samples.
    %
    %   OUTPUT
    %       noise    generated samples.
    %   some part of implementation found in http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
    %
    narginchk(2, 2);
    
    alpha = noiseDataSet.alpha;
    beta  = noiseDataSet.beta;
    
    if alpha == 1
        noise = -log(1 - rand(1, count)) * beta;
        return
    end
    
    alphaLessThenOne = alpha < 1;
    
    if alpha < 1
        alpha = alpha+1;
    end
    
    gamma = alpha-1;
    eta = sqrt(2.0 * alpha - 1.0);
    c = 0.5 - atan(gamma / eta) / pi;
    
    y(count) = 0;
    
    for k = 1:count,
        aux = -0.5;
        while aux < 0
            y(k) = -0.5;
            
            while y(k) <= 0
                u = rand(1, 1);
                y(k) = gamma + eta * tan(pi*(u - c) + c - 0.5);
            end
            
            v = -log(rand(1, 1));
            aux = v + log( 1.0 + ( (y(k)-gamma) / eta)^2 ) + gamma * log( y(k) / gamma ) - y(k) + gamma;
        end
    end
    
    if alphaLessThenOne
        noise = y .* beta .* (rand(1, count)).^(1.0 / (alpha - 1));
    else
        noise = y .* beta;
    end
end
