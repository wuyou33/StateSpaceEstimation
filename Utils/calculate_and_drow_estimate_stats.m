function calculate_and_drow_estimate_stats( errors, time_data, param_names, filter_type, model_name, units )
    % calculate_and_drow_estimate_stats.
    % Calculate mean and dquare-root standard deviation of filter estimation errors and plot results into separate figures.
    %
    %   void calculate_and_drow_estimate_stats( errors, time_data, legend, filter_type )
    %
    %   INPUT
    %       errors          matrix of error arrays: [MxPxN], where
    %                           M - number of runs to capture statistics;
    %                           P - number of estimated params;
    %                           N - number of samples (length of simaltion time);
    %       time_data       <TimeExt> struct which describe simulation time;
    %       param_names     names of estimated params;
    %       filter_type     filter name;
    %       model_name      inference data model name;
    %       units           measurement units of every param.
    %
    narginchk(6, 6);
    
    params_count = size(errors, 2);
    
    for i = 1:params_count
        param_name  = param_names(i);
        error_arr   = squeeze(errors(:, i, :));
        
        mean_x_err  = mean(error_arr);
        std_x_err   = std(error_arr);
        rmse_x_err = sqrt( mean( (error_arr).^2) );
        
        figure('Color', 'w');
        
        clf;
        p1 = plot(time_data.TimeInHour(2:end), mean_x_err, 'b', 'LineWidth', 1);
        hold on;
        p2 = plot(time_data.TimeInHour(2:end), mean_x_err - std_x_err, 'r', 'LineWidth', 1);
        hold on;
        p3 = plot(time_data.TimeInHour(2:end), mean_x_err + std_x_err, 'r', 'LineWidth', 1);
        hold off;
        grid on;
        
        legend([p1 p2 p3], char(strcat('mean rmse of', {' '}, param_name)), '- 1 sigma', '+ 1 sigma');
        xlabel('time, hours');
        ylabel(strcat(param_name, {' error, '}, units(i)));
        title({char(strcat(filter_type, ': Nonlinear Time Variant State Estimation of', {' '}, param_name, {' '})); ...
            char(strcat('(', model_name, ')'))});
        
        figure();
        clf;
        plot(time_data.TimeInHour(2:end), rmse_x_err, 'LineWidth', 1);
        grid on;
        xlabel('time, hours');
        ylabel(strcat({'RMSE '}, param_name, {', '}, units(i)));
        title(strcat(filter_type, ': RMSE of', {' '}, param_name, {' ('}, model_name, ')'));
        
        fprintf('RMSE %s = %d\n', char(param_name), rmse_x_err(end));
    end
end
