function skewMatrix = skew(vector)
    if 3 ~= size(vector,1),
        error('SCREWS:skew','vector must be 3x1')
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

    %   skewMatrix(1,1) = 0;
    %   skewMatrix(2,2) = 0;
    %   skewMatrix(3,3) = 0;
end