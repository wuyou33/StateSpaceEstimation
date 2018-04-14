function [ qout ] = quaternion_normalize( q )
    %  quaternion_normalize. Normalize a quaternion.
    %
    %   [ qout ] = quaternion_normalize( q )
    %
    %   Calculates the normalized quaternion, qout, for a given quaternion, q.  Input q is an 4-by-N matrix containing quaternions.
    %   qout returns an 4-by-N matrix of normalized quaternions. Each element of q must be a real number.
    %   Additionally, Q has its scalar number as the first column.
    %
    %   INPUT
    %       q   quaternion.
    %
    %   OUTPUT
    %       qout    normalized quaternion.
    %
    [m, n] = size(q);
    
    if m == 1
        error('q should be 4-by-N matrix or 4-by-1 vector')
    end
    
    for index = n:-1:1
        t = q(:, index);
        qnorm = norm(t, 2);
        qout(:, index) = t / qnorm;
    end
end
