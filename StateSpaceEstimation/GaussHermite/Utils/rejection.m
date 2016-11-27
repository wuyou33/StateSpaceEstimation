function [ reject ] = rejection( order, dim, indices )
    % rejection. Include or exclude current point set (by indices) in sparse Gauss-Hermite point set or not
    %
    %   [ reject ] = rejection( order, dim, indices )
    %
    %   INPUT
    %       order       order of Gauss-Hermite quadrature rule;
    %       dim         dimension of filtration (estimation) issue;
    %       indices     indeces.
    %
    %   OUTPUT
    %       reject  define indeces should be rejected in sparse Gauss-Hermite rule or not.
    %
    narginchk(3, 3);
    
    reject = sum(indices) <= order + dim - 1;
end
