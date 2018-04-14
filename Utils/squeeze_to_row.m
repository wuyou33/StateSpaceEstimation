function [ output_matrix ] = squeeze_to_row( input_matrix )
    % squeeze_to_row. Remove singleton dimensions from input matrix. If result is a one-dimensional array, then it would be a row vector.
    %
    %   [ output_matrix ] = squeeze_to_rowMatrix( input_matrix )
    %
    %   INPUT
    %       input_matrix   input matrix.
    %
    %   OUTPUT
    %       output_matrix   output matrix.
    %
    output_matrix = squeeze(input_matrix);
    dim = length(size(output_matrix));
    
    if dim < 3 && size(output_matrix, 2) == 1
        output_matrix = output_matrix';
    end
end
