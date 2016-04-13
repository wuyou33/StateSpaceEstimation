classdef WienerProcess < handle
    %WIENERPROCESS:
    %   A standard Wiener process (often called Brownian motion) on the interval $
    %   [0,T] is a random variable W(t) that depends continuously on $t in [0,T] and satisfies the following
    %   W(0) = 0
    %   0<= s <= t <= T
    %   dW = dt?.5 * N(0, 1)
    %   where N(0, 1) - gaussian noise
    
    %%
    properties (Access = private)
        mu;         % mean of process
        sigma;      % standard deviation
        dimension;  % dimension of process
    end
    
    %%
    methods (Access = public)
        function obj = WienerProcess(mu, sigma)
            if (length(mu) ~= length(sigma))
                error(' [WienerProcess] dimensions of mu and sigma does not agree ');
            end
            obj.dimension  = length(mu);
            obj.mu         = mu;
            obj.sigma      = sigma;
        end
        
        function simalationResults = simulate(this, sampleTime, sampleNumber)
            dw = sqrt(sampleTime) * repmat(this.sigma, 1, sampleNumber) .* randn(this.dimension, sampleNumber) + ...
                repmat(this.mu, 1, sampleNumber);   % increments
            
            simalationResults = cumsum(dw, 2)';
        end
    end
    
end
