function [ initialVal ] = loadInitialOrbit()
    % loadInitialOrbit. Load initial values of orbit.
    %
    %   OUTPUT
    %       [position; velocity]
    %       where:
    %           position - is a position in [km];
    %           velocity - is a velocity in [km / sec].
    %
    
    % GEO dreyf
    %     initialVal = [21588.3309009214; 42686.0207815723; 0; -2.57597028572808; 1.30278948238825; 0; 1; 0; 0; 0];
    
    % GEO2
    %     initialVal = [-21223.9926714100; -42868.3395565589; 0; 2.58697266398439; -1.28080278894642; 0; 1; 0; 0; 0];
    
    % GEO3
    %     initialVal = [-212234.9926714100; -4286.3395565589; 0; 2.58697266398439; -1.28080278894642; 0; 1; 0; 0; 0];
    
    % LO sat
    %     initialVal = [-6158.34755458333; -4063.45976435815; 0; 3.98653859590107; -6.04177022463943; 1.27633791914791; 1; 0; 0; 0];
    %
    
    % from gnss imitator
    % initialVal = [-6151.40172374421; -4073.96636606523; 2.22082127802645; 3.99717361592190; -6.03478370645393; 1.27632841244172; ...
    %         0.999997677569155; -0.000682877248575766; 0.00204467595732109; -1.35556652892711e-05];
    
    %%
    % todo: defined in ICRF-j2000, need to convert. ICRF-j2000 has a center as Solar System baricenter.
    % Hence require to move to the center of Coordinate syten to the center of Eath.
    
    %% Deep Impact Flyby/EPOXI Spacecraft (DI) (Start=2010-08-19, Stop=2010-09-18)
    % Deep Impact Flyby/EPOXI Spacecraft (DI). The EPOXI mission flyby reveals that comet Hartley 2's rocky ends spew out tons of
    % golf-ball to basketball-size fluffy ice particles, while the smooth middle area is more like what was observed on comet
    % Tempel 1 with water evaporating below the surface and percolating out through the dust.
    %
    epoxiVal = [1.263992563405297E+08; -7.120149330454805E+07; -6.710445880653515E+06; 1.561401205555831E+01; 2.767115733787828E+01; -1.009836672338032E+00];
    earthVal = [1.245753369794310E+08; -8.490915173563674E+07; 3.958015387110412E+03; 1.625451518298891E+01; 2.451202313450900E+01; -1.613990348628747E-03];
    initialVal = [epoxiVal(1:3) - earthVal(1:3); epoxiVal(4:6) - earthVal(4:6); 1; 0; 0; 0];
    
    %%  ExoMars TGO (Start=2016-07-22, Stop=2016-08-21)
    % ExoMars Mission objective is to search for evidence of methane and other trace
    % atmospheric gases that could be signatures of active biological or geological
    % processes, and to test technologies for future ESA contributions to Mars missions.
    % exoMarsVal = [2.944743436025662E+07;-1.932504435404160E+08;-1.260832133814412E+07;2.439494954830838E+01;4.676046384236798E-01;1.985076097147882E-01];
    % earthVal = [7.515763759387881E+07;-1.319970004198263E+08;-1.984076868406683E+04;2.545848216543012E+01;1.451279398144172E+01;4.671180059201419E-04];
    % initialVal = [exoMarsVal(1:3) - earthVal(1:3); exoMarsVal(4:6) - earthVal(4:6); 1; 0; 0; 0];
    
    %% NAVSTAR-68 (Start=2016-07-22, Stop=2016-08-21)
    %     navStarVal = [7.517195636922367E+07;-1.319756848263827E+08;-2.623423162838817E+04;2.343915251306972E+01;1.663284558568122E+01;2.559157449971319E+00];
    %     earthVal = [7.515763759387881E+07;-1.319970004198263E+08;-1.984076868406683E+04;2.545848216543012E+01;1.451279398144172E+01;4.671180059201419E-04];
    %     initialVal = [navStarVal(1:3) - earthVal(1:3); navStarVal(4:6) - earthVal(4:6); 1; 0; 0; 0];
    
end
