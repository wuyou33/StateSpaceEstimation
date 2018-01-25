function [ outIndex ] = multinomialResample( weights )
    % Multinomial resampling method for particle filters families.
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
    
    m = length(weights);
    q = cumsum(weights);
    q(m) = 1;
    
    outIndex = zeros(size(weights));
    i = 1;
    while (i <= m)
        sampl = rand(1, 1);  % (0,1]
        j = 1;
        
        while (q(j) < sampl)
            j = j+1;
        end;
        
        outIndex(i) = j;
        i = i+1;
    end
end
