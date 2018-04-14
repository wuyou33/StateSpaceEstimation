function pdf_evolution_plot_2(iterations, trueState, legend, legend_mu, timeData)
    figure();
    
    dim = size(iterations, 2);
    
    for i = 1 : dim
        for j = 1 : timeData.SimulationNumber
            pts = squeeze(iterations(:, i, j)) - trueState(:, i, j);
            [x_pdf, x_i] = ksdensity(pts, 'function', 'pdf', 'Kernel', 'box');
            t = timeData.TimeInHour(j) * ones(size(x_i));
            
            subplot(dim, 1, i);
            plot3(t, x_i, x_pdf, 'g', 'LineWidth', 0.5);
            hold on;
        end
        
        grid on;
        xlabel('time, hours');
        ylabel(strcat(legend(i), legend_mu(i)));
        zlabel('pdf');
        title(strcat('Evolution pdf of ', legend{i}));
    end
end
