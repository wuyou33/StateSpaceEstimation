function [ x, w ] = GaussHermite(n)
    %   GaussHermite
    %
    %   [ x, w ] = GaussHermite(n)
    %
    
    i   = 1 : n-1;
    a   = sqrt(i / 2);
    cm  = diag(a, 1) + diag(a, -1);
    
    
    [v, l]   = eig(cm);
    [x, ind] = sort(diag(l));
    
    v = v(:, ind)';
    w = sqrt(pi) * v(:, 1).^2;
    
    x = x' * sqrt(2);
    w = w' / sum(w);
    
end
