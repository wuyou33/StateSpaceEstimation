function [ points, weights ] = cubature_quadrature_points( stateDim, order )
    % cubature_quadrature_points. Draw cubature-quadrature points for cubature-quadrature Kalman filter.
    %
    %   Cubature points calculated as intersection of unit hyper-sphere and its axes.
    %   Quadrature points calculated as solution of Chebyshev-Laguerre polynoms with order n' and a = (n / 2 - 1).
    %
    %   [ points, weights ] = cubature_quadrature_points( stateDim, order )
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
    cubaturePoints   = intersect_unit_vector_hyper_sphere(stateDim);
    [ quadraturePoints, w ] = laguerre_quadrature_rule(order, alpha);
    
    num = 2*stateDim*order;
    weights = zeros(1, num);
    points = zeros(stateDim, num);
    
    for i = 1:order
        first = 2*stateDim*(i-1)+1;
        last  = 2*stateDim*i;
        
        points(:, first:last) = sqrt(2*quadraturePoints(i))*cubaturePoints;
        weights(first:last) = w(i)/(2*stateDim*gamma(stateDim*0.5));
    end
end
