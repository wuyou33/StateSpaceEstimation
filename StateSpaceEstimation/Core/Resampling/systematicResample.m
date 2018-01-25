function [ outIndex ] = systematicResample( weights )
    % Systematic resampling method for particle filters families.
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
    
    t = linspace(0, 1 - 1/n, n) + rand(1) / n;
    t(n+1) = 1;
    
    i=1;
    j=1;
    outIndex = zeros(size(weights));
    while (i <= n)
        if (t(i) < q(j))
            outIndex(i) = j;
            i = i+1;
        else
            j = j+1;
        end
    end
end
