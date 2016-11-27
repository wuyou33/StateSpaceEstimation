function [ output ] = mergeCells( firstArray, secondArray )
    % mergeCells. Merge two input one-dimensional cell arrays into one one-dimensional cell array.
    %
    %   [ output ] = mergeCells( firstArray, secondArray )
    %
    %   INPUT
    %       firstArray  	one-dimensional cell array;
    %       secondArray 	one-dimensional cell array.
    %
    %   OUTPUT
    %       output - one-dimensional cell array, which contains all elements from both input arrays (can have duplicates).
    %
    if ~iscell(firstArray)
        error('[ mergeCells::firstArray ] must be a one-dimensional cell array');
    end
    
    if ~iscell(secondArray)
        error('[ mergeCells:secondArray ] must be a one-dimensional cell array');
    end
    
    output = {firstArray{:}; secondArray{:}};
end
