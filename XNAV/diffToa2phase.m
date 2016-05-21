function [ phase ] = diffToa2phase( xRaySources, diffToa )
% diffToa2phase: Convert difference between toa to spaceshift and ssb to
% phase
% INPUT:
%       xRaySources - array of x-ray sources, every item should be instance of the XRaySource;
%       diffToa     - difference between toa to spaceshift and ssb (solar system baricenter);
% OUTPUT:
%       phase       - array of phase for every xRaySource;
%%     
    dimension = length(xRaySources);    
    invPeriods   = zeros(1, dimension);
    for i = 1:dimension        
        x = xRaySources(i);
        invPeriods(i) = 2*pi/x.period;
    end
    
    [capacity, ~] = size(diffToa);
    
    absPhase = diffToa.*repmat(invPeriods, capacity, 1);
    phase = absPhase;
%     phase = ang180(absPhase);
%     phase = unwrap(absPhase);
end