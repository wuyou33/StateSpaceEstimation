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
    %%
    narginchk(2, 2);
    
    %%
    particlesNum = length(weights);
    outIndex = zeros(1, particlesNum);
    
    weightsRes = particlesNum*weights;
    fix_kind = fix(weightsRes);
    
    % residual number of particles to sample
    numToSample = particlesNum - sum(fix_kind);
    
    if numToSample
        weightsRes = (weightsRes - fix_kind) / numToSample;
        cumDist = cumsum(weightsRes);
        
        u = fliplr(cumprod(rand(1, numToSample) .^ (1 ./ (numToSample:-1:1))));
        
        j = 1;
        for i = 1:numToSample
            while u(1, i) > cumDist(1, j)
                j = j + 1;
            end
            
            fix_kind(1, j) = fix_kind(1,j)+1;
        end
    end
    
    index = 1;
    
    for i = 1 : particlesNum
        if fix_kind(1, i) > 0
            for j = index : index + fix_kind(1, i) -1
                outIndex(j) = inIndex(i);
            end
        end
        
        index = index + fix_kind(1, i);
    end
end
