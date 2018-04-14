function [ centres, options, post ] = kmeans2( centres, data, options )
    % kmeans2. Algorithm to set the centres of a cluster model. The matrix data represents the data which is being clustered, with each row
    %	corresponding to a vector. The sum of squares error function is used. The point at which a local minimum is achieved is returned as
    %	centres.  The error value at that point is returned in options.squaredDistanceCluster.
    %
    %   Details:
    %       http://stats.stackexchange.com/questions/133656/how-to-understand-the-drawbacks-of-k-means
    %       https://en.wikipedia.org/wiki/K-means_clustering
    %       https://github.com/umath92/MachineLearning_5--k_means-kernel_k_means_gaussian_mixture_model/blob/master/kmeans.m
    %
    %	[ centres, options, post, errlog ] = kmeans( centres, data, options )
    %
    %   INPUT
    %       centres 	array of cluster centers;
    %       data    	matrix data represents the data which is being clustered, with each row corresponding to a vector;
    %       options 	algorithm options (details see bellow).
    %
    %   OUTPUT
    %       centres 	array of cluster centers;
    %       options 	options with changed (or not changed) squaredDistanceCluster value;
    %       post    	characterized points in cluster.
    %
    %	also returns the cluster number (in a one-of-N encoding) for each data point in 'post' and a log of the error values after each cycle in 'errlog'.
    %
    %   The optional parameters have the following interpretations.
    %
    %	options.logError                is set to 1 thenwarning messages are displayed. If this option is 0, then nothing is displayed.
    %
    %	options.absPrecision            is a measure of the absolute precision required for the value of CENTRES at the solution.
    %                                   If the absolute difference between the values of CENTRES between two successive steps is less than options.absPrecision,
    %                                   then this condition is satisfied.
    %
    %	options.precision               is a measure of the precision required of the error function at the solution.  If the absolute difference between
    %                                   the error functions between two successive steps is less than options.precision, then this condition is satisfied.
    %                                   Both this and the previous condition must be satisfied for termination.
    %
    %	options.iterationNumber         is the maximum number of iterations.
    %
    %   options.initializeFromData      if true then centres and posteriors need to be initialised from data
    %
    %   options.squaredDistanceCluster  is total squared distance from cluster centres
    %
    %   options.visualize               draw results.
    %
    narginchk(3, 3);
    
    logError            = options.logError;
    initializeFromData  = options.initializeFromData;
    precision           = options.precision;
    absPrecision        = options.absPrecision;
    iterationNumber     = options.iterationNumber;
    
    [ndata, dataDim] = size(data);
    [ncentres, dim]  = size(centres);
    
    if dataDim ~= dim;
        error('[ kmeans ] Data dimension does not match dimension of centres');
    end
    
    if ncentres > ndata
        error(' [ kmeans] More centres than data');
    end
    
    if initializeFromData == 1
        perm = randperm(ndata);
        perm = perm(1:ncentres);
        
        % assign first ncentres (permuted) data points as centres
        centres = data(perm, :);
    end
    
    id = eye(ncentres);
    
    for i = 1 : iterationNumber
        prevCentres = centres;
        
        % calculate posteriors based on existing centres
        dist = euclidean_distance(data, centres);
        
        % assign each point to nearest centres
        [minvals, index] = min(dist, [], 2);
        post = id(index, :);
        
        numPoints = sum(post, 1);
        % adjust the centres based on new posteriors
        for j = 1:ncentres
            if (numPoints(j) > 0)
                idx = (find(post(:, j)));
                centres(j, :) = sum(data(idx, :), 1) / numPoints(j);
            end
        end
        
        % error value is total squared distance from cluster centres
        err = sum(minvals);
        if logError > 0
            fprintf(1, 'Cycle %4d; Error %11.11f; \n', i, err);
        end
        
        % test for termination
        if i > 1 && (max( max( abs(centres - prevCentres) ) ) < absPrecision) && (abs(prevErr - err) < precision)
            options.squaredDistanceCluster = err;
            return;
        end
        
        prevErr = err;
    end
    
    % if we get here, then we haven't terminated in the given number of iterations.
    options.squaredDistanceCluster = err;
end
