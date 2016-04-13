function [ q ] = quaternionMultiply( q1, q2 )
%%  quaternionMultiply Calculate the product of two quaternions.
%   N = quaternionMultiply( Q, R ) calculates the quaternion product, N, for two
%   given quaternions, Q and R.  Inputs Q and R can be either 1-by-4 matrices 
%   containing a single 1-by-4 quaternion.  N returns an 
%   1-by-4 matrix of quaternion products. Each element of Q and R must be a
%   real number.  
%%
    % Calculate vector portion of quaternion product
    % vec = s1*v2 + s2*v1 + cross(v1,v2)
    
    vec = [q1(1).*q2(2) q1(1).*q2(3) q1(1).*q2(4)] + ...
            [q2(1).*q1(2) q2(1).*q1(3) q2(1).*q1(4)]+...
            [ q1(3).*q2(4)-q1(4).*q2(3) ...
            q1(4).*q2(2)-q1(2).*q2(4) ...
            q1(2).*q2(3)-q1(3).*q2(2)];

    % Calculate scalar portion of quaternion product
    % scalar = s1*s2 - dot(v1,v2)
    scalar = q1(1).*q2(1) - q1(2).*q2(2) - q1(3).*q2(3) - q1(4).*q2(4);
        
    q = [scalar  vec];   
end