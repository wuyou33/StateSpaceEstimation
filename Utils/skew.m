function [ skewMatrix ] = skew(vector)
    % skew. Calculated skew symmetric matrix from input vector.
    %
    %   [ skewMatrix ] = skew(vector)
    %
    %   INPUT
    %       vector  input 3D vector.
    %
    %   OUTPUT
    %       skewMatrix  corresponded skew symmetric matrix.
    %
    if 3 ~= size(vector,1),
        error('[ skew::skew ] vector must be 3x1')
    end
    
    if isnumeric(vector),
        skewMatrix = zeros(3,3);
    end
    
    skewMatrix(1,2) = -vector(3);
    skewMatrix(1,3) =  vector(2);
    skewMatrix(2,3) = -vector(1);
    
    skewMatrix(2,1) =  vector(3);
    skewMatrix(3,1) = -vector(2);
    skewMatrix(3,2) =  vector(1);
end
