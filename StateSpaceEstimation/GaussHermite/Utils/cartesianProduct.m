function [ z ] = cartesianProduct( varargin )
    % combinations - All combinations
    %   z = combinations(x1, x2, x3, ..., xN) returns all combinations of the elements
    %   in the arrays x1, x2, ..., and xN. z is P-by-N matrix is which P is the product
    %   of the number of elements of the N inputs. This functionality is also
    %   known as the Cartesian Product.
    %% error checking
    narginchk(1, Inf);
    if any(cellfun('isempty', varargin)); error('[ cartesianProduct:EmptyInput] One of more empty inputs result in an empty output.'); end
    %%
    args = varargin{1};
    num = length(args);
    ii = num : -1 : 1;
    
    [z{ii}] = ndgrid(args{ii});
    % concatenate
    z = reshape(cat(num+1, z{:}), [], num);
end
