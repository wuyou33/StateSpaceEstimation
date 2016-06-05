function [ initialVal ] = loadInitialOrbit()
%% Load initial values of orbit
% OUTPUT:
%   [position; velocity]
%   where: position - is a position in [km];
%          velocity - is a velocity in [km / sec].
%%
%     initialVal = [21588.3309009214; 42686.0207815723; 0; -2.57597028572808; 1.30278948238825; 0; 1; 0; 0; 0]; % GEO dreyf
    % initialVal = [-21223.9926714100; -42868.3395565589; 0; 2.58697266398439; -1.28080278894642; 0; 1; 0; 0; 0]; % HEO2
    initialVal = [-6158.34755458333; -4063.45976435815; 0; 3.98653859590107; -6.04177022463943; 1.27633791914791; 1; 0; 0; 0]; % LO sat
    % todo: several deep space trajectory
end