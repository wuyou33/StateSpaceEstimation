function [ output ] = merge_cells( first_array, second_array )
    % merge_cells. Merge two input one-dimensional cell arrays into one one-dimensional cell array.
    %
    %   [ output ] = merge_cells( first_array, second_array )
    %
    %   INPUT
    %       first_array  	one-dimensional cell array;
    %       second_array 	one-dimensional cell array.
    %
    %   OUTPUT
    %       output - one-dimensional cell array, which contains all elements from both input arrays (can have duplicates).
    %
    if ~iscell(first_array)
        error('[ merge_cells::first_array ] must be a one-dimensional cell array');
    end
    
    if ~iscell(second_array)
        error('[ merge_cells:second_array ] must be a one-dimensional cell array');
    end
    
    output = {first_array{:}; second_array{:}};
end
