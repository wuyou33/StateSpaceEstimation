function [ points, weigths ] = laguerreQuadratureRule( order, alpha )
    % laguerreQuadratureRule. Calculate cubature points using Gauss-Laguerre quadratures.
    %
    %   [ points, weigths ] = laguerreQuadratureRule( order, alpha )
    %
    %   Determines the abscisas (points) and weights (weights) for the Gauss-Laguerre quadrature of order n > 1, on the interval [0, +infinity].
    %   This is due to the fact that the companion matrix (of the n'th degree Laguerre polynomial) is now constructed as a symmetrical
    %   matrix, guaranteeing that all the eigenvalues (roots) will be real.
    %
    %   For detailed information please see:
    %       https://en.wikipedia.org/wiki/Classical_orthogonal_polynomials#Chebyshev_polynomials
    %
    %   Possible implementation via matlab function (more expensive and require additional computation of weights):
    %       syms z;
    %       y = vpa(solve(laguerreL(order, z), z));
    %
    %   Building the companion matrix CM.
    %   CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial under consideration.
    %   Moreover, CM will be constructed in such a way that it is symmetrical.
    %
    %   INPUT
    %       order   n'th degree Laguerre polynomial (subscript index of polynomial);
    %       alpha   parameter of Laguerre polynomial (superscript index polynomial), alpha = n / 2 - 1, where n - dimension of state space.
    %
    %   OUTPUT
    %       points      set of generated points;
    %       weights     set of corresponding weights for every point.
    %
    i   = 1:order;
    a   = (2*i-1) + alpha;
    b   = sqrt( i(1:order-1) .* ((1:order-1) + alpha) );
    cm  = diag(a) + diag(b, 1) + diag(b, -1);
    
    % Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
    
    [v, l]        = eig(cm);
    [points, ind] = sort(diag(l));
    
    v       = v(:, ind)';
    weigths = gamma(alpha+1) .* v(:, 1).^2;
end
