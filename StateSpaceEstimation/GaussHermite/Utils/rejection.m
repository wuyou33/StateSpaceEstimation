function reject = rejection(order, dim, indices)
    % rejection Include or exclude current point set (by indices) in sparse Gauss-Hermite point set or not
    reject = sum(indices) <= order + dim - 1;
end
