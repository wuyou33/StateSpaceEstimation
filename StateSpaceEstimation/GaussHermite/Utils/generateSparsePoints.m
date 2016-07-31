function [ points, weights ] = generateSparsePoints( order, indices, pointSI, manner)
    % generateSparsePoints. Generate subset of point for sparse Gauss-Hermite algorithm.
    %   [ points, weights ] = generateSparsePoints( order, indices, pointSI, manner)
    %
    %   INPUT
    %       order     - order of Gauss-Hermite polynomial;
    %       indices   - indeces, for which points and weights should be generated;
    %       pointSI   - set of point which determine that point must be generated or not;
    %       manner    - increase manner, L, 2L-1, or something else.
    %
    %   OUTPUT
    %       points  - generated points
    %       weights - generated weights for every point.
    %%
    dim = numel(indices);
    nq = sum(indices) - dim;
    
    if nq >= order-dim && nq <= order-1
        [points, weights] = generateSingleDimPoint(indices(1), pointSI(1), manner);
        
        for i = 2 : dim
            [pt, w] = generateSingleDimPoint(indices(i), pointSI(i), manner);
            
            points  = extend(points, pt, numel(w), size(weights, 2));
            weights = extend(weights, w, numel(w), size(weights, 2));
        end
        
        if size(weights, 1) ~= 1
            weights = prod(weights);
            weights = weights.*( factorial(dim-1) / ( factorial(order-1-nq) * factorial((dim-1)-(order-1-nq)) ) ) * (-1)^(order-1-nq);
        end
    else
        points = [];
        weights = [];
    end
end

function extended = extend( orig, additonal, repOrig, repAdd )
    orig = repmat(orig, 1, repOrig);
    ext = repmat(additonal, repAdd, 1);
    extended = [orig; ext(:)'];
end