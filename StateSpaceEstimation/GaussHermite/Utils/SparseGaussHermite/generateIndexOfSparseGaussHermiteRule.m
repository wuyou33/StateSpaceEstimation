function [ res ] = generateIndexOfSparseGaussHermiteRule(accuracyLevel, dimension)
    %   generateIndexOfSparseGaussHermiteRule
    %
    %   [ res ] = generateIndexOfSparseGaussHermiteRule(accuracyLevel, dimension)
    %
    nindex   = generateNextIndexOfSparseGaussHermiteRule(accuracyLevel, ones(dimension, 1), []);
    nnIndex  = [ones(dimension,1), nindex];
    tmpIndex = [];
    
    while 1
        for i = 1 : size(nindex, 2)
            tmpIndex = generateNextIndexOfSparseGaussHermiteRule(accuracyLevel, nindex(:, i), nnIndex);
            nnIndex  = [nnIndex, tmpIndex];
        end
        
        if isempty(tmpIndex)==1
            break;
        end
        
        nindex = nnIndex;
    end
    
    res = nnIndex;
end
