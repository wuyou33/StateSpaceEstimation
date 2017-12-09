function [ newParticles ] = multivariateSmoothResampling( particles, weights )
    % multivariateSmoothResampling. Multivariate smooth resampling for SIR (sequential importance sampling).
    %
    %   [ new_particles ] = multivariateSmoothResampling( particles, weights )
    %
    %   Performs the resampling stage of the SIR algorithm in order (number of samples) steps.
    %
    %   INPUT
    %       particles   particles at step k;
    %       weights 	corresponding weights at step k.
    %
    %   OUTPUT
    %       newParticles    resampled particles.
    %
    %%
    narginchk(2, 2);
    
    %%
    particlesNum = length(weights);
    stateDim = size(particles, 2);
    
    [P, D] = eig(particles'*(bsxfun(@times, 1 / particlesNum, particles)));
    
    D = diag(D);
    
    vectors = bsxfun(@times, P, sqrt(D)');
    
    orthogonalizedParticles = bsxfun(@rdivide, particles*vectors, D');
    newParticles = zeros(particlesNum, stateDim);
    
    for j = 1:stateDim
        idxOut = sortrows([orthogonalizedParticles(:,j) weights], 1);
        newParticles(:, j) = univariateSmoothResampling(idxOut(:, 2), idxOut(:, 1), particlesNum);
    end
    
    newParticles = newParticles*(vectors');
    newParticles = newParticles';
end
