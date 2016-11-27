function [ flag ] = isInteger( x )
    % isInteger. Check that input data is integer number.
    %
    %   [ flag ] = isInteger( x )
    %
    %   INPUT
    %       x   input arg.
    %
    %   OUTPUT
    %       flag    result of checking.
    %
    if isempty(x)
        error('Input cannot be empty');
    elseif ~isnumeric(x)
        error('Input must be numeric');
    end
    
    flag = isfinite(x) && ( x == floor(x) );
end
