function [ sources ] = loadXRaySources(count)
    % loadXRaySources. Load x-ray sources for X-Nav
    %
    %   [ sources ] = loadXRaySources(count)
    %
    %   INPUT
    %       count 	number of required x-ray sources. should be positive integer between 1 and 7.
    %
    %    OUTPUT
    %       sources     loaded X-Ray sources.
    %
    if ~isnumeric(count)
        error('[ loadXRaySources ] count should be positive integer between 1 and 8');
    end
    
    if (count < 1 || count > 7)
        error('[ loadXRaySources ] count should be positive integer between 1 and 8');
    end
    
    p1 = XRaySource('B0531+21', ... name it's Crub Pulsar
        33.51e-3,... period
        1.535,... intensity
        75,... raError
        60,... decError
        7.7e-3,... gSource
        3.3e-3,... gBackgr
        184.6,... galacticLon
        -5.78,... galacticLat
        2.0 ... distance
        );
    
    %     p2 = XRaySource('B1821-24', ... name
    %         3.05e-3,... period
    %         1.93e-4,... intensity
    %         0.9,... raError
    %         12,... decError
    %         1.9e-3,... gSource
    %         1.1e-3,... gBackgr
    %         7.8,... galacticLon
    %         -5.58,... galacticLat
    %         4.9 ... distance
    %     );
    
    %     p3 = XRaySource('B1937+21', ... name
    %         1.56e-3,... period
    %         4.99e-5,... intensity
    %         0.012,... raError
    %         0.14,... decError
    %         7.41e-4,... gSource
    %         3.29e-4,... gBackgr
    %         57.5,... galacticLon
    %         -0.29,... galacticLat
    %         5.0 ... distance
    %     );
    
    %     p4 = XRaySource('J0218+4232', ... name
    %         2.32e-3,... period
    %         6.65e-5,... intensity
    %         150,... raError
    %         100,... decError
    %         1.18e-4,... gSource
    %         1.18e-4,... gBackgr
    %         139.5,... galacticLon
    %         -17.53,... galacticLat
    %         5.8 ... distance
    %     );
    
    p5 = XRaySource('B0540-69', ... name
        50.4e-3,... period
        5.15e-3,... intensity
        550,... raError
        500,... decError
        2.55e-3,... gSource
        2.55e-3,... gBackgr
        279.7,... galacticLon
        -31.5,... galacticLat
        0.8 ... distance
        );
    
    p6 = XRaySource('J1814-338', ... name
        3.18e-3,... period
        9.97e-2,... intensity
        1,... raError -- really unknown
        1,... decError -- really unknown
        1.27e-2,... gSource
        1.27e-2,... gBackgr
        358.75,... galacticLon
        -7.59,... galacticLat
        1.5 ... distance -- approximatelly
        );
    
    p7 = XRaySource('J1808-369', ... name
        2.49e-3,... period
        3.29e-1,... intensity
        1,... raError -- really unknown
        1,... decError -- really unknown
        1.26e-4,... gSource
        1.268e-4,... gBackgr
        355.39,... galacticLon
        -8.15,... galacticLat
        0.5 ... distance -- approximatelly
        );
    
    %     allSources = [p1; p2; p3; p4; p5; p6; p7];
    allSources = [p1; p5; p6; p7];
    
    if count > length(allSources)
        n = length(allSources);
    else
        n = count;
    end
    
    sources = allSources(1:n);
end
