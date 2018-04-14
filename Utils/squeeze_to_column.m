function [ output_matrix ] = squeeze_to_column( input_matrix )
    % squeeze_to_column Remove singleton dimensions from input matrix. If result is a one-dimensional array, then it would be a column vector.
    %
    % 	[ output_matrix ] = squeezeToRowMatrix( input_matrix )
    %
    %   INPUT
    %       input_matrix  input matrix.
    %
    %   OUTPUT
    %       output_matrix   output matrix.
    %
    output_matrix = squeeze(input_matrix);
    dim = length(size(output_matrix));
    
    if dim < 3 && size(output_matrix, 1) == 1
        output_matrix = output_matrix';
    end
end
