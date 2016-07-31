function [ idx ] = generateNextIndices( order, indices, prevIndices )
    % generateNextIndices
    %   [ idx ] = generateNextIndices( order, indices, prevIndices )
    %
    %   INPUT:
    %       order       - order of Gauss-Hermite polynomial;
    %       indices     - current indices;
    %       prevIndices - previous generated indices.
    %
    %   OUTPUT:
    %       idx - generated indeces.
    %%
    dim = numel(indices);
    idx = zeros(dim,dim);
    
    for i = 1 : dim
        idx(:, i) = indices;
        idx(i, i) = idx(i, i) + 1;
        
        if rejection(order, dim, idx(:,i)) ~= 1
            idx(:, i) = 0;
        end
    end
    
    %% erase empty column;
    idx(:, sum(abs(idx), 1) == 0) = [];
    
    for i = 1 : size(prevIndices, 2)
        residual = sum( abs( idx - cvecrep(prevIndices(:, i), size(idx, 2)) ) );
        idx(:, residual == 0) = 0;
    end
    
    idx(:, sum(abs(idx), 1) == 0) = [];
end
