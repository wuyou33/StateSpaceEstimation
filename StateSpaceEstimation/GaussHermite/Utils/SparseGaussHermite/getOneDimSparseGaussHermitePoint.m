function [ gridPoints, gridWeights ] = getOneDimSparseGaussHermitePoint(index, pointSet, manner)
    %   getOneDimSparseGaussHermitePoint
    %
    %   [ gridPoints, gridWeights ] = getOneDimSparseGaussHermitePoint(index, pointSet, manner)
    %
    if pointSet == 1
        
        if index == 0
            gridPoints=[];
            gridWeights=[];
        else
            switch (manner)
                case 1
                    n = index;
                case 2
                    n = 2*index-1;
                case 3
                    n = 2^index-1;
                otherwise
                    error(' [getOneDPoint]:manner not supported ');
            end
            
            [ gridPoints, gridWeights ] = GaussHermite(n);
        end
    end
    
end
