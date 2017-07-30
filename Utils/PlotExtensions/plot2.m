function plot2( x, y, titleText, legendText, yLabelText, xLabelText )
    % plot2. Plots figure.
    %
    %   Pplots vector Y versus vector X. If X or Y is a matrix, then the vector is plotted versus the rows or columns of the matrix, whichever line up.
    %   If X is a scalar and Y is a vector, disconnected line objects are created and plotted as discrete points vertically at X.
    %
    %   INPUT
    %       x             vector;
    %       y             vector or matrix;
    %       titleText     title;
    %       legendText    legend;
    %       yLabelText    y axis label;
    %       xLabelText    x axis label.
    %
    narginchk(5, 6);
    
    if nargin < 6
        xLabelText = 'time, sec';
    end
    
    if isvector(y)
        plot(x, y, 'LineWidth', 1);
        grid on;
        hold on;
        title(titleText);
        legend(legendText);
        xlabel(xLabelText);
        ylabel(yLabelText);
    elseif ismatrix(y)
        for i = 1:size(y, 1)
            plot(x, y(i, :), 'LineWidth', 1);
            grid on;
            hold on;
        end
        title(titleText);
        legend(legendText);
        xlabel(xLabelText);
        ylabel(yLabelText);
    else
        error(' [plot2::y] must be vector or matrix');
    end
end
