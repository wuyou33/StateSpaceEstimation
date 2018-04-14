function [ out_index ] = residual_resample( input_index, weights )
    % residual_resample. Residual resampling for SIR (sequential importance sampling).
    %
    %   http://stackoverflow.com/questions/14091662/bootstrapping-hierarchical-multilevel-data-resampling-clusters/14131934#14131934
    %   http://robotics.stackexchange.com/questions/479/particle-filters-how-to-do-resampling
    %
    %   [ out_index ] = residual_resample( inIndex, weights )
    %
    %   Performs the resampling stage of the SIR algorithm in order (number of samples) steps.
    %
    %   INPUT
    %       input_index     (r-vector) input particle indices;
    %       weights         (r-vector) normalised importance weights.
    %
    %   OUTPUT
    %       out_index  - resampled indices.
    %
    %
    narginchk(2, 2);
    
    particles_count = length(weights);
    
    out_index = zeros(1,particles_count);
    
    % first integer part
    weights_residual = particles_count * weights;
    integer_weights_kind = fix(weights_residual);
    
    residual_particles_count = particles_count-sum(integer_weights_kind);
    
    if residual_particles_count        
        weights_residual = (weights_residual - integer_weights_kind) / residual_particles_count;
        cum_dist = cumsum(weights_residual);
        
        % generate N (N = residual_particles_count) ordered random variables uniformly distributed in [0, 1]
        u = fliplr(cumprod(rand(1, residual_particles_count) .^ (1 ./ (residual_particles_count : -1 :1))));
        j = 1;
        
        for i = 1 : residual_particles_count
            while (u(1, i) > cum_dist(1, j))
                j = j+1;
            end
            integer_weights_kind(1, j) = integer_weights_kind(1, j) + 1;
        end        
    end
    
    index = 1;
    for i = 1:particles_count
        if (integer_weights_kind(1, i) > 0)
            upper_index = index + integer_weights_kind(1, i) - 1;
            for j = index : upper_index
                out_index(j) = input_index(i);
            end;
        end;
        
        index = index + integer_weights_kind(1, i);
    end
end
