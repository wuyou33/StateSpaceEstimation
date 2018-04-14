function [ outIndex ] = resample( resample_method, weights, particles_count )
    % Resample procedure Summary for particle filters families.
    % Supports following methods
    %   - residual;
    %   - multinomial;
    %   - stratified;
    %   - systematic;
    %   - residual2 (alternative residual implemetation which supports kitagawa sampling method).
    %
    %   [ outIndex ] = resample( resample_method, weights, particles_count )
    %
    %   INPUT:
    %       resample_method     resample method which will be used to resample particles at step k;
    %       weights             arrray of weights of particle filter at step k;
    %       particles_count     total count of particels.
    %
    %   OUTPUT:
    %       outIndex    resampled indices.
    %
    narginchk(3, 3);
    
    switch (resample_method)
        case 'residual'
            outIndex  = residual_resample(1:particles_count, weights);
        case 'multinomial'
            outIndex = multinomial_resample(weights);
        case 'stratified'
            outIndex = stratified_resample(weights);
        case 'systematic'
            outIndex = systematic_resample(weights);
        case 'residual2'
            outIndex = residual_resample_2(weights, rand(size(weights)));
        case 'residual-kitagawa'
            outIndex = residual_resample_2(weights, rand(1));
        otherwise
            error('[ resample ] unknown inference_model.resampleMethod type.');
    end
end
