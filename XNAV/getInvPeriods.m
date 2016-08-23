function [ invPeriods ] = getInvPeriods( xRaySources )
    % GETINVPERIODS Calculate following expression 2*pi / Period for every xRaySource.
    %
    %   [ invPeriods ] = getInvPeriods( xRaySources )
    %
    %   INPUT:
    %       xRaySources - array of the X-Ray sources, every item should be instance of the XRaySource.
    %
    %   OUTPUT:
    %       invPeriods  - result of 2*pi / Period.
    
    dimension = length(xRaySources);    
    invPeriods   = zeros(1, dimension);
    for i = 1:dimension
        x = xRaySources(i);
        invPeriods(i) = x.TwoPiOnPeriod;
    end
end
