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
            outIndex  = residualResample(1:particles_count, weights);
        case 'multinomial'
            outIndex = multinomialResample(weights);
        case 'stratified'
            outIndex = stratifiedResample(weights);
        case 'systematic'
            outIndex = systematicResample(weights);
        case 'residual2'
            outIndex = residualResample2(weights, rand(size(weights)));
        case 'residual-kitagawa'
            outIndex = residualResample2(weights, rand(1));
        otherwise
            error('[ resample ] unknown inference_model.resampleMethod type.');
    end
end
