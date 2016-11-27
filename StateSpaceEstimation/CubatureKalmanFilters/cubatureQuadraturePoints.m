function [ points, weights ] = cubatureQuadraturePoints( stateDim, order )
    % cubatureQuadraturePoints. Draw cubature-quadrature points for cubature-quadrature Kalman filter.
    %
    %   Cubature points calculated as intersection of unit hyper-sphere and its axes.
    %   Quadrature points calculated as solution of Chebyshev-Laguerre polynoms with order n' and a = (n / 2 - 1).
    %
    %   [ points, weights ] = cubatureQuadraturePoints( stateDim, order )
    %
    %   INPUT
    %       stateDim    state space dimenstion;
    %       order       order of the laguerre quadrature rule.
    %
    %   OUTPUT
    %       points      matrix of cubature points;
    %       weights     array of corresponded weights.
    %
    narginchk(2, 2);
    
    alpha    = stateDim*0.5 - 1;
    cubaturePoints   = intersectUnitVectorHyperSphere(stateDim);
    [ quadraturePoints, weights ] = laguerreQuadratureRule(order, alpha);
    
    weights = weights/(2*stateDim*gamma(stateDim*0.5));
    
    points = zeros(stateDim, 2*stateDim*order);
    for i = 1:order
        points(:, 2*stateDim*(i-1)+1 : 2*stateDim*i) = sqrt(2)*cubaturePoints*quadraturePoints(i);
    end
end
