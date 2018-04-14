function m = column_vector_replicate(v, c)
    % column_vector_replicate. Column vector replicate.
    %
    %    m = column_vector_replicate(v, c)
    %
    %   Replicates a Nx1 dimensional column vector v, c times to generate a NxC dimensional matrix m.
    %
    %   INPUT
    %       v   vector;
    %       c   count of replication.
    %
    %   OUTPUT
    %       m   matrix of replicated vector.
    %
    if isempty(v)
        m = zeros(0, c);
    else
        m = v(:, ones(c, 1));
    end
end
