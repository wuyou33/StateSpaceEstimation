function draw_error_on_iteration( time_data, errors_1, errors_2, filter_type, labels, legends_1, legends_2, model_name )
    % draw_error_on_iteration
    subplot(2, 1, 1);
    plot(time_data.TimeInHour(2:end), errors_1, 'LineWidth', 1);
    xlabel('time, hours');
    ylabel(labels{1});
    legend(legends_1);
    hold on;
    grid on;
    title({char(strcat(filter_type, {': Nonlinear State Estimation errors'})); char(strcat('(', model_name, ')'))});
    
    subplot(2, 1, 2);
    plot(time_data.TimeInHour(2:end), errors_2, 'LineWidth', 1);
    xlabel('time, hours');
    ylabel(labels{2});
    legend(legends_2);
    hold on;
    grid on;
    drawnow
end
