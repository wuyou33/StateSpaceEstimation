function plotClaster( x, k )
    % plotClaster. Plot 2d/3d samples of different classes with different colors.
    %
    %   plotClaster( x, clusterCnt )
    %
    %   INPUT
    %       x   data (matrix); k   clusters count.
    %
    %%
    %# K-means clustering
    %# (K: number of clusters, G: assigned groups, C: cluster centers)
    [G, C] = kmeans(x, k, 'distance', 'sqEuclidean', 'start', 'sample');
    
    %# show points and clusters (color-coded)
    clr = lines(k);
    figure();
    hold on;
    scatter3(x(:, 1), x(:, 2), x(:, 3), 36, clr(G,:), 'Marker',  '.');
    scatter3(C(:, 1), C(:, 2), C(:, 3), 100, clr, 'Marker', 'o', 'LineWidth', 3);
    
    hold off;
    view(3),
    axis vis3d,
    box on,
    rotate3d on
    
    xlabel('x'),
    ylabel('y'),
    zlabel('z');
    
    axis equal;
    grid on;
    hold off;
end
