function plot2( x, y, titleText, legendText, yLabelText, xLabelText )
% PLOT2 plots vector Y versus vector X. If X or Y is a matrix,
%   then the vector is plotted versus the rows or columns of the matrix,
%   whichever line up.  If X is a scalar and Y is a vector, disconnected
%   line objects are created and plotted as discrete points vertically at X.
% INPUT:
%       x           - vector
%       y           - vector
%       titleText   - title
%       legendText  - legend,
%       yLabelText  - y axis label
%       xLabelText  - x axis label
    
    if nargin < 6
        xLabelText = 'time, sec';
    end
    
    plot(x, y, 'LineWidth', 1);
    grid on;
    hold on;
    title(titleText);
    legend(legendText);
    xlabel(xLabelText);
    ylabel(yLabelText);
end
