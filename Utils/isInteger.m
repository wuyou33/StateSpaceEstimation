function [ flag ] = isInteger( x )
    % ISINTEGER Check that input data is integer number.
    if isempty(x)
        error('Input cannot be empty');
    elseif ~isnumeric(x)
        error('Input must be numeric');
    end
    
    flag = isfinite(x) && ( x == floor(x) );
end
