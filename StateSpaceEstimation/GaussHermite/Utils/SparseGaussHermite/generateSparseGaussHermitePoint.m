function [ points, weights ] = generateSparseGaussHermitePoint(accuracyLevel, currentIndex, pointSet, manner)
    %   generateSparseGaussHermitePoint
    %   [ points, weights ] = generateSparseGaussHermitePoint(accuracyLevel, currentIndex, pointSet, manner)
    %
    
    dim = numel(currentIndex);
    q = sum(currentIndex) - dim;
    
    if q >= accuracyLevel - dim && q <= accuracyLevel - 1
        
        [points, weights] = getOneDimSparseGaussHermitePoint(currentIndex(1),pointSet(1),manner);
        
        for i = 2 : dim
            [npt, nw] = getOneDimSparseGaussHermitePoint(currentIndex(i), pointSet(i), manner);
            
            num_npt = numel(nw);
            num_pt  = size(weights,2);
            
            points = repmat(points, 1, num_npt);
            pt_add = repmat(npt, num_pt, 1);
            pt_add = (pt_add(:))';
            points = [points; pt_add];
            
            weights = repmat(weights,1,num_npt);
            w_add   = repmat(nw, num_pt,1);
            w_add   = (w_add(:))';
            weights = [weights;w_add];
        end
        
        if size(weights, 1)
            weights = prod(weights);
            weights = weights.*( ...
                factorial(dim - 1) / ( factorial(accuracyLevel-1-q) * factorial((dim-1) - (accuracyLevel - 1 -q)) ) ...
                ) * (-1)^(accuracyLevel - 1 -q);
        end
        
    else
        points=[];
        weights=[];
    end
    
end