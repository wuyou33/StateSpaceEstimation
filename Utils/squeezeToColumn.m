function [ outputMatrix ] = squeezeToColumn( inputMatrix )
    % squeezeToColumn Remove singleton dimensions from input matrix. If result is a one-dimensional array, then it would be a column vector.
    %
    % 	[ outputMatrix ] = squeezeToRowMatrix( inputMatrix )
    %
    %   INPUT
    %       inputMatrix  input matrix.
    %
    %   OUTPUT
    %       outputMatrix   output matrix.
    %
    outputMatrix = squeeze(inputMatrix);
    dim = length(size(outputMatrix));
    
    if dim < 3 && size(outputMatrix, 1) == 1
        outputMatrix = outputMatrix';
    end
end
