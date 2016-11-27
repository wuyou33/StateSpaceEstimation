function [ matrix ] = elrep( v, m, n )
    % elrep. Replicate every element of input vector x-by-y times.
    %
    %   Where x = m is number of replication by row and y = n / dimension(v) is a number of replication by column
    %   Detailed explanation goes here
    %
    %   [ matrix ] = elrep( v, m, n )
    %
    %   INPUT
    %       v 	vector Nx1 dimension;
    %       m 	number of row replication;
    %       n 	number of column replication.
    %
    %   OUPUT
    %       matrix 	matrix dimension (m x N)-by-n.
    %
    if isempty(v)
        matrix = zeros(m, n);
    else
        matrix = zeros(m, n);
        x = n / length(v);
        for i = 1:length(v)
            matrix(1:m, x*(i-1)+1 : x*i) = repmat(v(i), m, x);
        end
    end
end
