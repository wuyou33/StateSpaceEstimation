function [ points, weights ] = generateSingleDimPoint( order, pointSI, manner )
    % generateSingleDimPoint.
    %
    %   [ points, weights ] = generateSingleDimPoint( order, pointSI, manner )
    %
    %   INPUT
    %       order     - order of Gauss-Hermite polynomial;
    %       pointSI   - indicate that point should be generated or not;
    %       manner    - increase manner, L, 2L-1, or something else.
    %
    %   OUTPUT
    %       points  - generated points
    %       weights - generated weights for every point.
    %
    narginchk(3, 3);
    
    if pointSI ~= 1 || order == 0
        points = [];
        weights=[];
        return;
    end
    
    switch manner
        case 1
            [points, weights] = gaussHermiteRule(order);
        case 2
            [points, weights] = gaussHermiteRule(2*order-1);
        case 3
            [points, weights] = gaussHermiteRule(2^order-1);
        otherwise
            error('[ generateSingleDimPoint::manner ] manner not supported. manner should be 1, 2 or 3');
    end
end
