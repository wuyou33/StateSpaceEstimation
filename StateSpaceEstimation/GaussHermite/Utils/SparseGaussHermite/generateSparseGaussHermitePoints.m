function [ points, weights ] = generateSparseGaussHermitePoints(accuracyLevel, dimension, manner)
    % generateSparseGaussHermitePoints
    %
    %   [ points, weights ] = generateSparseGaussHermitePoints(accuracyLevel, dimension, manner)
    %
    %   INPUT:
    %       accuracyLevel: accuracy level;
    %       dimension: state spase dimension;
    %       manner: increase manner: L, 2*L-1, 2^L-1.
    %
    %   OUTPUT:
    %       points      - array of generated points;
    %       weights     - vecoter of generated weights.
    %
    
    indexSet = generateIndexOfSparseGaussHermiteRule(accuracyLevel, dimension);
    
    points  = [];
    weights = [];
    pointSI = ones(1, dimension);
    
    for i = 1 : size(indexSet, 2)
        [tmpPoints, tmpWeights] = generateSparseGaussHermitePoint(accuracyLevel, indexSet(:, i), pointSI, manner);
        
        if isempty(points)
            points  = [points,  tmpPoints];
            weights = [weights, tmpWeights];
        else
            indexToDel = [];
            npt  = numel(weights);
            
            for j = 1 : numel(tmpWeights)
                residual = abs(points - cvecrep(tmpPoints(:, j), npt));
                fi = find(sum(residual) < 1e-6);
                
                if ~isempty(fi)
                    weights(fi)  = tmpWeights(j) + weights(fi);
                    indexToDel   = [indexToDel, j];
                end
            end
            
            tmpPoints(:, indexToDel)  = [];
            tmpWeights(indexToDel)    = [];
            
            points  = [points,  tmpPoints];
            weights = [weights, tmpWeights];
        end
    end
    
end
