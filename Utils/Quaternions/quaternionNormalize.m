function [ qout ] = quaternionNormalize( q )
    %%  QUATNORMALIZE Normalize a quaternion.
    %   N = quaternionNormilize( Q ) calculates the normalized quaternion, N, for a
    %   given quaternion, Q.  Input Q is an 4-by-N matrix containing
    %   quaternions.  N returns an 4-by-N matrix of normalized quaternions.
    %   Each element of Q must be a real number.  Additionally, Q has its
    %   scalar number as the first column.
    %%
    [m, n] = size(q);
    
    if m == 1; error('q should be 4-by-N matrix or 4-by-1 vector'); end
    
    for index = n:-1:1
        t = q(:, index);
        qnorm = norm(t, 2);
        qout(:, index) = t / qnorm;
    end
    
%     qMod = norm(q, 2);
%     qout2 = q / qMod;
%     
%     if ~isequal(qout, qout2)
%         error('something went wrong');
%     end
end
