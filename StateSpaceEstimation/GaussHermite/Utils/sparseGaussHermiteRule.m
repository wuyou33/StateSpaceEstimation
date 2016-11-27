function [ points, weights ] = sparseGaussHermiteRule( order, dimension, manner )
    % sparseGaussHermiteRule. Draw points and corresponded weights for sparse Gauss-Hermite rule.
    %
    %   Generate points and weights according to Sparse Gauss-Hermite quadrature rule.
    %   Weight function is chosen to be the standard Gaussian density with zero mean and unit variance N(0; I).
    %   The interval of interest is chosen to be (-infinity; +infinity).
    %   According to the fundamental theorem of Gauss-Hermite quadrature, the quadrature points are chosen to be the zeros of the m-th order Hermite polynomial.
    %   Since the zeros of the Hermite polynomials are distinct,
    %   it is noteworthy that the determinant of the coefficient matrix in is the well known Vandermonde's determinant that is nonzero.
    %   For an m-point quadrature scheme, the resulting quadrature rule is exact for all polynomials of degree  <= 2m - 1
    %
    %   [ points, weights ] = sparseGaussHermiteRule( order, dimension, manner )
    %
    %   INPUT
    %       order     	order of Gauss-Hermite polynomial;
    %       dimension 	dimension;
    %       manner    	increase manner, L, 2L-1, or something else.
    %
    %   OUTPUT
    %       points  	generated points
    %       weights 	generated weights for every point.
    %
    narginchk(3, 3);
    
    indices = generateIndices(order, dimension);
    points  = [];
    weights = [];
    pointSI = ones(1, dimension);
    
    for i = 1:size(indices, 2)
        [pt, w] = generateSparsePoints(order, indices(:, i), pointSI, manner);
        
        numPoints = numel(weights);
        if isempty(points) == 1
            points = [points pt];
            weights = [weights w];
        else
            idx = [];
            
            for j = 1 : numel(w)
                residual = abs(points-repmat(pt(:,j),1,numPoints));
                fi = find(sum(residual)<1e-6);
                
                if numel(fi) > 1
                    error('[ sparseGaussHermiteRule ] Something went wrong');
                end
                
                if isempty(fi) ~= 1
                    weights(fi) = w(j) + weights(fi);
                    idx = [idx j];
                end
            end
            
            pt(:, idx) = [];
            w(idx) = [];
            
            points = [points pt];
            weights = [weights w];
        end
    end
end
