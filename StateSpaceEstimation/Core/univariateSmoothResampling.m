function [ newParticles ] = univariateSmoothResampling( weights, particles, newParticlesNumber )
    % multivariateSmoothResampling. Multivariate smooth resampling for SIR (sequential importance sampling).
    %
    %   [ new_particles ] = multivariateSmoothResampling( particles, weights )
    %
    %   Performs the resampling stage of the SIR algorithm in order (number of samples) steps.
    %
    %   INPUT
    %       particles              particles at step k;
    %       weights                corresponding weights at step k;
    %       newParticlesNumber     number of new particles to generate.
    %
    %   OUTPUT
    %       newParticles    resampled particles.
    %
    %   based on https://github.com/DynareTeam/particles/blob/master/src/univariate_smooth_resampling.m
    %%
    narginchk(3, 3);
    
    %%
    cnt = length(particles);
    
    lambdaTilde = [ (.5*weights(1)); (.5*(weights(1:cnt-1)+weights(2:cnt))); (.5*weights(cnt)) ];    
    lambdaBar = cumsum(lambdaTilde);
    
    u = rand(1, 1) ;
    
    newParticles = zeros(newParticlesNumber,1) ;
    
    rj = 0;
    i = 1;
    j = 1;
    
    while i <= newParticlesNumber
        u_j = ( i-1 + u)/newParticlesNumber;
        
        while u_j > lambdaBar(j)
            rj = j;
            j = j+1;
        end
        
        if rj == 0
            newParticles(i) = particles(1);
        elseif rj==cnt
            newParticles(i) = particles(cnt);
        else
            u_star = (u_j - lambdaBar(rj))./lambdaTilde(rj+1);
            newParticles(i) = (particles(rj+1) - particles(rj))*u_star + particles(rj);
        end
        
        if isnan(newParticles(i))
            display('a');
        end
        i = i+1;
    end
end
