function [ phase ] = diffToa2phase( invPeriods, diffToa )
    % diffToa2phase. Convert difference between toa to spaceshift and ssb to phase.
    %
    %   [ phase ] = diffToa2phase( invPeriods, diffToa )
    %
    %   INPUT
    %       invPeriods    array of 2*pi derived on period of x-ray source;
    %       diffToa       difference between toa to spaceshift and ssb (solar system baricenter).
    %
    %   OUTPUT
    %       phase   array of phase for every xRaySource.
    %
    capacity = size(diffToa, 2);
    
    absPhase = diffToa.*cvecrep(invPeriods', capacity);
    %     phase = absPhase;
    phase = ang180(absPhase);
    %     phase = unwrap(absPhase);
end
