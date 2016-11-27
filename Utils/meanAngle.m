function [ angleMean ] = meanAngle(angleList, limits, dim)
    % meanAngle.  Calculates the mean angle, correcting for angle discontinuities
    %
    %   [ angleMean ] = meanAngle(angleList, limits, dim)
    %
    %   Calculates the mean angle of the input vector/matrix angleList taking into account discontinuities.
    %   If dim is specified, the mean over dimension DIM is taken. angleList can be up to 2D.
    %
    %   INPUT
    %       angleList   input vector or matrix of angles (in degrees), up to two dimensions;
    %       limits      definition of angles: [-pi pi] or [0 2*pi] degrees;
    %       dim         dimension over which to operate.
    %
    %   OUTPUT
    %       angleMean 	mean of angleList (in radians).
    %
    %   Pascal de Theije, v2.0, 13 April 2005
    %   Copyright (c) 2005, TNO-D&V
    %   All Rights Reserved
    %
    
    % Check if input angles match with specified angle limits.
    limitsInDegree = 2*pi*limits;
    
    % if ~isempty(find((angleList < min(limitsInDegree))|(angleList > max(limitsInDegree))))
    %   disp('meanangle.m: Input angles do not match with limits.')
    %   return
    % end
    
    % If the dimension is not specified, operate over the first dimension.
    if (nargin == 2)
        dim = 1;
    end
    
    %
    % Transpose matrix, if necessary, so that all operations have to be done
    % over the second dimension.
    %
    if (dim == 1)
        angleList = angleList.';
    end
    
    %
    % In order to find the best estimate of the mean angle, calculate the
    % in-product of an arbitrary vector [a b] and all specified angles, and find
    % that direction that gives the highest inproduct. The in-product is given by
    % [a b].*[sum(cos(ANGLES)) sum(sin(ANGLES))]. The derivative of this w.r.t.
    % 'a' is set to zero, which gives the solution a=sqrt(C^2/(C^2+S^2)).
    %
    C  = sum(cos(angleList*pi/180),2);
    S  = sum(sin(angleList*pi/180),2);
    CS = C.^2 + S.^2;
    % Calculate vector that gives highest inproduct.
    a  = sqrt(C.^2./CS);
    b  = sqrt(S.^2./CS);
    
    %
    % From the way 'a' and 'b' are solved, they can be either positive or negative,
    % while in practice only one of these combinations is the true one. Find the
    % combination of 'a' and 'b' that gives the lowest in-product. This angle
    % is used to 'split' the angle circle.
    %
    temp = (C.*a)*[-1 1 -1 1] + (S.*b)*[-1 -1 1 1];
    [~, ind] = max(temp,[],2);
    ind2 = find(ind==1 | ind==3);
    a(ind2) = -a(ind2);
    ind2 = find(ind==1 | ind==2);
    b(ind2) = -b(ind2);
    
    %
    % Find angle that should be used to 'split' the angle circle.
    %
    cut_angle = atan2(b,a) * 180 / pi;
    
    %
    % Split the angle circle at 'cut_angle'.
    %
    angleList = ang180(angleList - cut_angle*ones(1,size(angleList,2)));
    % Then take 'normal' mean.
    angleMean = mean(angleList,2);
    % Undo splitting angle circle.
    angleMean = ang180(angleMean + cut_angle);
    
    %
    % Transform output to angles between 0 and 360 degrees, if necessary.
    %
    if (limitsInDegree == [0 360])
        ind = find(angleMean < 0);
        angleMean(ind) = angleMean(ind) + 360;
    end
    
    %
    % Transpose matrix, if necessary.
    %
    if (dim == 1)
        angleMean = angleMean.';
    end
    
    angleMean = angleMean/2/pi;
end
