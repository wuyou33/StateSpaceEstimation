function [ m ] = rvecrep(v, c)
% RVECREP  Row vector replicate
%   M = rvecrep(V, C) Replicates a 1xN dimensional row vector V, C times to generate a CxN dimensional matrix M.
%%
    if isempty(v)
        m = zeros(c, 0);
    else
        m = v(ones(1, c), :);
    end
end
