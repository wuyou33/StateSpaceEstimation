function m = cvecrep(v, c)
%%
%    Column vector replicate
%    M = cvecrep(V, C) Replicates a Nx1 dimensional column vector V, C times to generate a NxC dimensional matrix M.
%%
    if isempty(v)
        m = zeros(0,c);
    else
        m = v(:, ones(c,1));
    end
end
