classdef ExpCorrGaussianProcess < handle
    % ExpCorrGaussianProcess. Exponentially correlated Gaussian random process.
    %   An algorithm is described which generates a sequence of random numbers r1, r2, . . . with the following two
    %   properties: (i) each individual ri is a Gaussian deviate with zero mean and unit variance; (ii) the autocorrelation
    %   function of the sequence decays exponentially with a predetermined decay time ? . A correlated random walk is
    
    properties (Access = private)
        corrTime;
        mean;
        sigma;
    end
    
    methods (Access = public)
        function obj = ExpCorrGaussianProcess(corrTime, mean, sigma)
            if ~isequal(size(corrTime), size(mean), size(sigma))
                error('[ ExpCorrGaussianProcess ] dimension mismatch');
            end
            
            obj.corrTime    = corrTime;
            obj.mean        = mean;
            obj.sigma       = sigma;
        end
        
        function res = simulate(this, sampleNumber)
            dim = length(this.mean);
            g = cvecrep(this.mean, sampleNumber) + cvecrep(this.sigma, sampleNumber) .* randn(dim, sampleNumber);
            
            res = zeros(dim, sampleNumber);
            res(:, 1) = g(:, 1);
            
            f = exp(-ones(dim, 1) ./ this.corrTime);
            ff = (ones(dim, 1) - f.^2).^0.5;
            
            for i = 2:sampleNumber
                res(:, i) = f.*res(:, i-1) + ff.*g(:, i);
            end
        end
    end
    
end
