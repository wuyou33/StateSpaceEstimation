function [ points ] = intersectUnitVectorHyperSphere( n )
    % intersectUnitVectorHyperSphere. Return intersection points between a unit hyper sphere and its axis.
    %   For additional details please see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    %
    %   [ points ] = intersectUnitVectorHyperSphere( n )
    %
    %   INPUT
    %       n - space dimension.
    %
    %   OUTPUT
    %       points - intersection points between a unit hyper sphere and its axis
    % 
    % unit sphere. for 3D: [x_center y_center z_center  R]
    sphere = [zeros(n), ones(n, 1)];    
    
    % unit vectors. for 3D: [x0 y0 z0  dx dy dz]
    line = [zeros(n), eye(n)];
        
    points = intersectLineHyperSphere(line, sphere, n);
end
