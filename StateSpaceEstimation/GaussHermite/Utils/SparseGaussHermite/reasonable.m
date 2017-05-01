function [ z ] = reasonable(accuracyLevel, dimension, index)
    %   reasonable
    %
    %   [ z ] = reasonable(accuracyLevel, dimension, index)
    %
    
    ci = sum(index);
    
    if ci <= accuracyLevel + dimension - 1
        z = 1;
    else
        z = 0;
    end
end
