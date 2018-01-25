function [ outIndex ] = residualResample2( weights, noise )
    % residualResample. Residual resampling for SIR (sequential importance sampling).
    %
    %   http://stackoverflow.com/questions/14091662/bootstrapping-hierarchical-multilevel-data-resampling-clusters/14131934#14131934
    %   http://robotics.stackexchange.com/questions/479/particle-filters-how-to-do-resampling
    %
    %   [ outIndex ] = residualResample2( particles, weights, noise )
    %
    %   Performs the resampling stage of the SIR algorithm in order (number of samples) steps.
    %
    %   INPUT
    %       weights       normalised importance weights;
    %       noise         noise matrix.
    %
    %   OUTPUT
    %       outIndex      resampled indices.
    %
    narginchk(2, 2);
    
    particlesNum = length(weights);
    
    switch length(noise)
        case 1
            kitagawaResampling = 1;
        case particlesNum
            kitagawaResampling = 0;
        otherwise
            error(['residualResample2: Unknown method! The size of the second argument (' inputname(3) ') is wrong.'])
    end
    
    jndx = 1 : particlesNum;
    outIndex = zeros(1, particlesNum);
    
    mWeights = particlesNum * weights;
    iWeights = fix(mWeights);
    
    % Compute the number of resample
    trialsNum = particlesNum - sum(iWeights);
    
    if trialsNum
        mWeights = (mWeights - iWeights)/trialsNum;
        empiricalCDF = cumsum(mWeights);
        
        if kitagawaResampling
            u = ( transpose(1:trialsNum) - 1 + noise(:) ) / trialsNum;
        else
            u = fliplr( cumprod(noise(1:trialsNum) .^ (1./(trialsNum:-1:1))) );
        end
        
        j = 1;
        for i = 1:trialsNum
            while (u(i) > empiricalCDF(j))
                j = j+1;
            end
            
            iWeights(j) = iWeights(j) + 1;
            if kitagawaResampling == 0
                j = 1;
            end
        end
    end
    
    k = 1;
    for i = 1:particlesNum
        if iWeights(i) > 0
            for j = k : k + iWeights(i) - 1
                outIndex(j) = jndx(i);
            end
        end
        
        k = k + iWeights(i);
    end
end
