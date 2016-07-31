function [ points, weights ] = cubatureQuadraturePoints( stateDim, order )
    % cubatureQuadraturePoints Summary of this function goes here
    %   Detailed explanation goes here
    %%
    alpha    = stateDim*0.5 - 1;
    cubaturePoints   = intersectUnitVectorHyperSphere(stateDim);
    [ quadraturePoints, weights ] = laguerreQuadratureRule(order, alpha);
    
    weights = weights/(2*stateDim*gamma(stateDim*0.5));
    
    points = zeros(stateDim, 2*stateDim*order);        
    for i = 1:order
        points(:, 2*stateDim*(i-1)+1 : 2*stateDim*i) = sqrt(2)*cubaturePoints*quadraturePoints(i);        
    end
end

