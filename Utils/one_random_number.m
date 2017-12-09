function [ num ] = one_random_number( lower, upper, exclude )
    % ONE_RANDOM_NUMBER. Generate random numbers between two numbers without any specific number
    %
    %   INPUT
    %       lower       lower bound;
    %       upper       upper bound;
    %       exclude     value to exclude.
    %
    %   OUTPUT
    %       num     generated random number;
    %
    narginchk(2, 3);
    
    if (nargin == 2)
        exclude = NaN;
    end
    
    num = exclude;
    
    while isnan(num) || any( num == exclude )
        num = lower + rand(1, 1) * (upper - lower);
    end
end
