function [ z ] = cartesianProduct( varargin )
    % cartesianProduct. Find all combinations for input set (varargin)
    %
    %   [ z ] = cartesianProduct( varargin )
    %
    %   z = combinations(x1, x2, x3, ..., xN) returns all combinations of the elements in the arrays x1, x2, ..., and xN.
    %   z is P-by-N matrix is which P is the product of the number of elements of the N inputs.
    %   This functionality is also known as the Cartesian Product.
    %
    %   INPUT
    %       varargin    input set (array of cells).
    %
    %   OUTPUT
    %       z   all combinations of the elements of input arrays.
    %
    narginchk(1, Inf);
    
    if any(cellfun('isempty', varargin))
        error('[ cartesianProduct:EmptyInput] One of more empty inputs result in an empty output.');
    end
    
    args = varargin{1};
    num = length(args);
    ii = num : -1 : 1;
    
    [z{ii}] = ndgrid(args{ii});
    
    % concatenate
    z = reshape(cat(num+1, z{:}), [], num);
end
