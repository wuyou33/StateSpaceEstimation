function [ indices ] = generateIndices( order, dimension )
    % generateIndices. Generate indices for sparse Gauss-Hermite rule.
    %   [ indices ] = generateIndices( order, dimension )
    %
    %   INPUT:
    %       order     - order of Gauss-Hermite polynomial;
    %       dimension - dimension;
    %   OUTPUT:
    %       indices - generated indeces.
    %%
    
    nextIndex = generateNextIndices(order, ones(dimension, 1), []);
    nnIndex = [ones(dimension, 1) nextIndex];
    tmpIndex = [];
    
    while 1
        for i = 1 : size(nextIndex, 2)
            tmpIndex = generateNextIndices(order, nextIndex(:, i), nnIndex);
            nnIndex = [nnIndex tmpIndex];
        end
        
        if isempty(tmpIndex) == 1
            break;
        end
    
        nextIndex = nnIndex;
    end
    
    indices = nnIndex;
end
