function squareRootMatrix = svdDecomposition(matrix)
    %% Calculate square root factor of matrix through Singular value decomposition
    
    %%
    [u, s, v] =  svd(matrix);
    squareRootMatrix = 0.5*(u + v) * sqrt(s);
end

