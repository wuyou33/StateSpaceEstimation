function [ gridPoints, gridWeights ] = get_one_dim_sparse_gauss_hermite_point(index, pointSet, manner)
    %   get_one_dim_sparse_gauss_hermite_point
    %
    %   [ gridPoints, gridWeights ] = get_one_dim_sparse_gauss_hermite_point(index, pointSet, manner)
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
            
            [ gridPoints, gridWeights ] = gauss_hermite(n);
        end
    end
    
end
