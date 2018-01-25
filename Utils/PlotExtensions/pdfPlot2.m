function pdfPlot2(iterations, trueState, legend)
    figure();
    
    dim = size(iterations, 2);
    
    for i = 1 : dim
        pts = squeeze(iterations(:, i, end)) - trueState(:, i, end);
        [x_pdf, x_i] = ksdensity(pts, 'function', 'pdf', 'Kernel', 'box');
        [x_pdf_norm, ~] = ksdensity(pts, 'function', 'pdf', 'Kernel', 'normal');
        
        subplot(dim, 1, i);
        plot(x_i, x_pdf, '--r', 'LineWidth', 2);
        hold on;
        
        subplot(dim, 1, i);
        plot(x_i, x_pdf_norm, '--b', 'LineWidth', 1);
        hold on;
        
        subplot(dim, 1, i);
        h1 = histogram(pts, 'Normalization', 'pdf');
        h1(1).FaceColor = [.8 .8 1];
        hold on;
        
        grid on;
        
        title(strcat(legend{i}, ' pdf'));
    end
end
