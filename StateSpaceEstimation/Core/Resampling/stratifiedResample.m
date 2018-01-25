function [ outIndex ] = stratifiedResample( weights )
    % Stratified resampling method for particle filters families.
    %
    %   [ outIndex ] = multinomialResample( weights )
    %
    %   INPUT:
    %       weights     arrray of weights of particle filter at step k.
    %
    %   OUTPUT:
    %       outIndex    resampled indices.
    %
    narginchk(1, 1);
    
    n = length(weights);
    q = cumsum(weights);
    
    for i = 1:n
        t(i) = rand(1, 1) / n + (i - 1) / n;
    end
    
    t(n+1) = 1;
    
    i = 1;
    j = 1;
    outIndex = zeros(size(weights));
    while ( i <= n )
        if (t(i) < q(j))
            outIndex(i) = j;
            i = i+1;
        else
            j = j+1;
        end
    end
end
