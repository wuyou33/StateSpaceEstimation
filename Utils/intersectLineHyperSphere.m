function [ points ] = intersectLineHyperSphere( line, sphere, n )
    % intersectUnitVectorHyperSphere Return intersection points between a line in hyper space and a hyper sphere.
    %   For additional details please see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    %
    %   [ points ] = intersectLineHyperSphere( line, sphere, n )
    %   INPUT
    %       line    - array of line coordinates and in hyper space, for instance [x0 y0 z0 ... n0 dx dy dz ... dn];
    %       sphere  - array of center points and radius of hyper sphere, for instance [x_center y_center z_center ... n_center  R];
    %       n       - dimension of hyper space.
    %       
    %   OUTPUT
    %       points - intersection points between a unit hyper sphere and its axis.
    %%
    tolerance = 1e-14;
    
    diffCenters = bsxfun(@minus, line(:, 1:n), sphere(:, 1:end-1));
    
    % equation coefficients
    a = sum(line(:, n:end) .* line(:, n:end), 2);    
    b = 2 * sum(bsxfun(@times, diffCenters, line(:, n+1:end)), 2);    
    c = sum(diffCenters.*diffCenters, 2) - sphere(:, end).*sphere(:, end);
    
    % solve equation
    delta = b.*b - 4*a.*c;
    
    % initialize empty results
    points = NaN * ones(2 * size(delta, 1), n);
    
    % proces couples with two intersection points
    inds = find(delta > tolerance);
    if ~isempty(inds)
        % delta positive: find two roots of second order equation
        u1 = (-b(inds) -sqrt(delta(inds))) / 2 ./ a(inds);
        u2 = (-b(inds) +sqrt(delta(inds))) / 2 ./ a(inds);
        
        % convert into 3D coordinate
        points(inds, :) = line(inds, 1:n) + bsxfun(@times, u1, line(inds, n+1:end));
        points(inds+length(delta), :) = line(inds, 1:n) + bsxfun(@times, u2, line(inds, n+1:end));
    end
    
    % process couples with one intersection point
    inds = find(abs(delta) < tolerance);
    if ~isempty(inds)
        % delta around zero: find unique root, and convert to 3D coord.
        u = -b(inds) / 2 ./ a(inds);
        
        % convert into 3D coordinate
        pts = line(inds, 1:end) + bsxfun(@times, u, line(inds, n+1:end));
        points(inds, :) = pts;
        points(inds+length(delta), :) = pts;
    end
end
