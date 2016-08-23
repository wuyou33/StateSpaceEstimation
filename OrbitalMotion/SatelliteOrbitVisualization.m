function [] = SatelliteOrbitVisualization(satellitePhaseState)

if (~isa(satellitePhaseState, 'SatellitePhaseSpace'))
    error('satellitePhaseState should be instance of the SatellitePhaseSpace');
end

grs80 = referenceEllipsoid('grs80','km');

figure('Renderer', 'opengl')

ax = axesm('globe', 'Geoid', grs80,...
    'Grid', 'off', ...
    'GLineWidth', 1,...
    'GLineStyle', '-'...
);

ax.Position = [0 0 1 1];
ax.XTickMode = 'auto';
ax.YTickMode = 'auto';
ax.ZTickMode = 'auto';
ax.XColor = [0 0.447058826684952 0.74117648601532];
ax.YColor = [0 0.447058826684952 0.74117648601532];
ax.ZColor = [0 0.447058826684952 0.74117648601532];
ax.Color = [0.941176474094391 0.941176474094391 0.941176474094391];
ax.PlotBoxAspectRatio = [360 320 120];

box(ax,'on');
grid(ax,'on');
hold(ax,'on');

view(3)

load topo
geoshow(topo, topolegend, 'DisplayType', 'texturemap')
demcmap(topo)

plot3(satellitePhaseState.TrajectoryX, ...
    satellitePhaseState.TrajectoryY,...
    satellitePhaseState.TrajectoryZ, ...
    '-black', ...
    'LineWidth', ...
    2 ...
);

figure;
hold on
plot3(satellitePhaseState.TrajectoryX, ...
    satellitePhaseState.TrajectoryY,...
    satellitePhaseState.TrajectoryZ, ...
    '-black', ...
    'LineWidth', 2 ...
);
legend('satellite trajectory');
grid on;
hold off; 