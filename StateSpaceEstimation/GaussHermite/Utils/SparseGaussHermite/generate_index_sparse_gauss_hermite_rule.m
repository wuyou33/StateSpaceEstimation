function [ res ] = generate_index_sparse_gauss_hermite_rule(accuracyLevel, dimension)
    %   genera_ieInd_sfSpar_geGau_hsHermi_reRule
    %
    %   [ res ] = generate_index_sparse_gauss_hermite_rule(accuracyLevel, dimension)
    %
    nindex   = generate_next_index_sparse_gauss_hermite_rule(accuracyLevel, ones(dimension, 1), []);
    nnIndex  = [ones(dimension,1), nindex];
    tmpIndex = [];
    
    while 1
        for i = 1 : size(nindex, 2)
            tmpIndex = generate_next_index_sparse_gauss_hermite_rule(accuracyLevel, nindex(:, i), nnIndex);
            nnIndex  = [nnIndex, tmpIndex];
        end
        
        if isempty(tmpIndex)==1
            break;
        end
        
        nindex = nnIndex;
    end
    
    res = nnIndex;
end
