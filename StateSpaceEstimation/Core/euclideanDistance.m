function [ distanceMatrix ] = euclideanDistance( left, right )
    % euclideanDistance. Takes two matrices of vectors and calculates the squared Euclidean distance between them.
    % Both matrices must be of the same column dimension. If Left is M-by-N matrix, and Right is L-by-N matrix, then the result is M-by-L matrix.
    % The I, Jth entry is the squared distance from the Ith row of Left to the Jth row of Right.
    %
    %   More about algorithm in
    %       http://stackoverflow.com/questions/23911670/efficiently-compute-pairwise-squared-euclidean-distance-in-matlab#
    %       http://math.stackexchange.com/questions/1236465/euclidean-distance-and-dot-product
    %
    %   [ distanceMatrix ] = euclideanDistance( left, right )
    %
    %   INPUT
    %       left    M-by-N matrix;
    %       right   M-by-L matrix.
    %
    %   OUPTUT
    %       distanceMatrix   M-by-L matrix where I, Jth entry is the squared distance from the Ith row of Left to the Jth row of Right.
    %
    narginchk(2, 2);
    
    [dimLeft1, dimLeft2] = size(left);
    [dimRight1, dimRight2] = size(right);
    
    if dimLeft2 ~= dimRight2
        error('[ euclideanDistance ] dimension mismatch');
    end
    
    distanceMatrix = (ones(dimRight1, 1) * sum(left.^2, 2)')' + ones(dimLeft1, 1) * sum(right.^2, 2)' - 2*( left*(right') );
    
    % reset all negative values to zero
    distanceMatrix(distanceMatrix < 0) = 0;
end
