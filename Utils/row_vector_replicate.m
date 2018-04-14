function [ m ] = row_vector_replicate(v, c)
    % row_vector_replicate. Row vector replicate.
    %   Replicates a 1xN dimensional row vector v, c times to generate a CxN dimensional matrix m.
    %
    %   [ m ] = row_vector_replicate(v, c)
    %
    %   INPUT
    %       v   input vector;
    %       c   number of replications.
    %
    %   OUTPUT
    %       m   matrix of replicated vector.
    %
    if isempty(v)
        m = zeros(c, 0);
    else
        m = v(ones(1, c), :);
    end
end
