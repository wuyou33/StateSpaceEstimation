function [ angle_mean ] = mean_angle(angle_list, limits, dim)
    % mean_angle.  Calculates the mean angle, correcting for angle discontinuities
    %
    %   [ angle_mean ] = mean_angle(angle_list, limits, dim)
    %
    %   Calculates the mean angle of the input vector/matrix angle_list taking into account discontinuities.
    %   If dim is specified, the mean over dimension DIM is taken. angle_list can be up to 2D.
    %
    %   INPUT
    %       angle_list   input vector or matrix of angles (in degrees), up to two dimensions;
    %       limits      definition of angles: [-pi pi] or [0 2*pi] degrees;
    %       dim         dimension over which to operate.
    %
    %   OUTPUT
    %       angle_mean 	mean of angle_list (in radians).
    %
    %   Pascal de Theije, v2.0, 13 April 2005
    %   Copyright (c) 2005, TNO-D&V
    %   All Rights Reserved
    %
    
    % Check if input angles match with specified angle limits.
    limitsInDegree = 2*pi*limits;
    
    % if ~isempty(find((angle_list < min(limitsInDegree))|(angle_list > max(limitsInDegree))))
    %   disp('mean_angle.m: Input angles do not match with limits.')
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
        angle_list = angle_list.';
    end
    
    %
    % In order to find the best estimate of the mean angle, calculate the
    % in-product of an arbitrary vector [a b] and all specified angles, and find
    % that direction that gives the highest inproduct. The in-product is given by
    % [a b].*[sum(cos(ANGLES)) sum(sin(ANGLES))]. The derivative of this w.r.t.
    % 'a' is set to zero, which gives the solution a=sqrt(C^2/(C^2+S^2)).
    %
    C  = sum(cos(angle_list*pi/180),2);
    S  = sum(sin(angle_list*pi/180),2);
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
    angle_list = unwrap_angle_from_min_pi_to_pi(angle_list - cut_angle*ones(1,size(angle_list,2)));
    % Then take 'normal' mean.
    angle_mean = mean(angle_list,2);
    % Undo splitting angle circle.
    angle_mean = unwrap_angle_from_min_pi_to_pi(angle_mean + cut_angle);
    
    %
    % Transform output to angles between 0 and 360 degrees, if necessary.
    %
    if (limitsInDegree == [0 360])
        ind = find(angle_mean < 0);
        angle_mean(ind) = angle_mean(ind) + 360;
    end
    
    %
    % Transpose matrix, if necessary.
    %
    if (dim == 1)
        angle_mean = angle_mean.';
    end
    
    angle_mean = angle_mean/2/pi;
end
