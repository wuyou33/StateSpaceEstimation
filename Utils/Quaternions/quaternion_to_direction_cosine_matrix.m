function [ dcm ] = quaternion_to_direction_cosine_matrix(q)
    %  quaternion_to_direction_cosine_matrix. Convert quaternion to direction cosine matrix.
    %
    %   [ dcm ] = quaternion_to_direction_cosine_matrix( q )
    %
    %   Calculates the direction cosine matrix, dcm, for a given quaternion, q. Input q is an 1-by-4 matrix containing M quaternions.
    %   N returns a 3-by-3 matrix of direction cosine matrices.  The direction cosine matrix performs the coordinate transformation of a vector in inertial
    %   axes to a vector in body axes. Each element of Q must be a real number.  Additionally, q has its scalar number as the first column.
    %
    %   INPUT
    %       q   quaternion.
    %
    %   OUTPUT
    %       dcm     direction cosine matrix.
    %
    qin = quaternion_normalize(q);
    
    dcm = zeros(3, 3);
    
    dcm(1,1) = qin(1).^2 + qin(2).^2 - qin(3).^2 - qin(4).^2;
    dcm(1,2) = 2.*(qin(2).*qin(3) + qin(1).*qin(4));
    dcm(1,3) = 2.*(qin(2).*qin(4) - qin(1).*qin(3));
    dcm(2,1) = 2.*(qin(2).*qin(3) - qin(1).*qin(4));
    dcm(2,2) = qin(1).^2 - qin(2).^2 + qin(3).^2 - qin(4).^2;
    dcm(2,3) = 2.*(qin(3).*qin(4) + qin(1).*qin(2));
    dcm(3,1) = 2.*(qin(2).*qin(4) + qin(1).*qin(3));
    dcm(3,2) = 2.*(qin(3).*qin(4) - qin(1).*qin(2));
    dcm(3,3) = qin(1).^2 - qin(2).^2 - qin(3).^2 + qin(4).^2;
    
end
