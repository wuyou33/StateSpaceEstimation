function loglog2( x, y, titleText, legendText, yLabelText, xLabelText )
    % loglog2. Plots figure in log scale.
    %
    %   Plots vector Y versus vector X. If X or Y is a matrix, then the vector is plotted versus the rows or columns of the matrix whichever line up.
    %   If X is a scalar and Y is a vector, disconnected line objects are created and plotted as discrete points vertically at X.
    %   Logarithmic (10) scale is used for both axis (the X-axis and the Y-axis).
    %
    %   INPUT
    %       x             vector;
    %       y             vector;
    %       titleText     title;
    %       legendText    legend;
    %       yLabelText    y axis label;
    %       xLabelText    x axis label.
    %
    narginchk(5, 6);
    
    if nargin < 6
        xLabelText = 'time, sec';
    end
    
    loglog(x, abs(y));
    grid on;
    hold on;
    title(titleText);
    legend(legendText);
    xlabel(xLabelText);
    ylabel(yLabelText);
end
