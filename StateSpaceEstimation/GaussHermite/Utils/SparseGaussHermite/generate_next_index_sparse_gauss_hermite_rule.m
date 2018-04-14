function [ res ] = generate_next_index_sparse_gauss_hermite_rule(accuracyLevel, index, oldIndex)
    %   generate_next_index_sparse_gauss_hermite_rule
    %
    %   generate_next_index_sparse_gauss_hermite_rule(accuracyLevel, index, oldIndex)
    %
    dim = numel(index);
    
    res = zeros(dim, dim);
    
    for i = 1 : dim
        res(:, i) = index;
        res(i, i) = res(i, i) + 1;
        
        if reasonable(accuracyLevel, dim, res(:,i)) ~= 1
            res(:, i) = 0;
        end
    end
    
    % erase empty column;
    res(:, sum(abs(res), 1) == 0) = [];
    nr   = size(res, 2);
    nold = size(oldIndex,2);
    
    for i = 1 : nold
        residual = abs(res-column_vector_replicate(oldIndex(:, i), nr));
        res(:, sum(residual) == 0) = 0;
    end
    
    res(:, sum(abs(res), 1) == 0) = [];
    
end
