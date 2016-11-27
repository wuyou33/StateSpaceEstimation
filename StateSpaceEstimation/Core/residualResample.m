function [ outIndex ] = residualResample( inIndex, weights )
    % residualResample. Residual resampling for SIR (sequential importance sampling).
    %
    %   http://stackoverflow.com/questions/14091662/bootstrapping-hierarchical-multilevel-data-resampling-clusters/14131934#14131934
    %   http://robotics.stackexchange.com/questions/479/particle-filters-how-to-do-resampling
    %
    %   [ outIndex ] = residualResample( inIndex, weights )
    %
    %   Performs the resampling stage of the SIR algorithm in order (number of samples) steps.
    %
    %   INPUT
    %       inIndex     (r-vector) input particle indices;
    %       weights 	(r-vector) normalised importance weights.
    %
    %   OUTPUT
    %       outIndex  - resampled indices.
    %
    narginchk(2, 2);
    
    partilcesNum = length(weights);
    outIndex = zeros(1, partilcesNum);
    
    % residual resampling
    weightsRes = partilcesNum * weights;
    nKind = fix(weightsRes);
    
    % residual number of particles to sample
    nRes = partilcesNum - sum(nKind);
    
    if nRes
        cumDist = cumsum((weightsRes - nKind) / nRes);
        
        % generate n result ordered random variables uniformly distributed in [0, 1]
        u = fliplr( cumprod(rand(1, nRes).^(1./(nRes : -1 : 1))) );
        
        j = 1;
        for i = 1:nRes
            while u(1, i) > cumDist(1, j)
                j = j + 1;
            end
            
            nKind(1, j) = nKind(1, j) + 1;
        end
    end
    
    % copy resampled trajectories
    index = 1;
    
    for i = 1 : partilcesNum
        if nKind(1, i) > 0
            endIndex = index + nKind(1, i) - 1;
            for j = index : endIndex
                outIndex(j) = inIndex(i);
            end
        end
        
        index = index + nKind(1, i);
    end
    
end
