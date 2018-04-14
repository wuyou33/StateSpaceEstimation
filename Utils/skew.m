function [ skew_matrix ] = skew(vector)
    % skew. Calculated skew symmetric matrix from input vector.
    %
    %   [ skew_matrix ] = skew(vector)
    %
    %   INPUT
    %       vector  input 3D vector.
    %
    %   OUTPUT
    %       skew_matrix  corresponded skew symmetric matrix.
    %
    if 3 ~= size(vector,1),
        error('[ skew::vector ] must be 3x1')
    end
    
    if isnumeric(vector),
        skew_matrix = zeros(3,3);
    end
    
    skew_matrix(1,2) = -vector(3);
    skew_matrix(1,3) =  vector(2);
    skew_matrix(2,3) = -vector(1);
    
    skew_matrix(2,1) =  vector(3);
    skew_matrix(3,1) = -vector(2);
    skew_matrix(3,2) =  vector(1);
end
