clear;
close all;
clc;

tic;
M = 3;
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWEO = Archive;
cost_MOGWEO = reshape([Archive_MOGWEO.Cost], M, [])';
load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWO = Archive;
cost_MOGWO = reshape([Archive_MOGWO.Cost], M, [])';
load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOEO = Archive;
cost_MOEO = reshape([Archive_MOEO.Cost], M, [])';

clear Archive;

if M == 2
    figure;
    ax = gca;
    ax.FontName = 'Times New Roman';
    plot(cost_MOGWO(:, 1), cost_MOGWO(:, 2), '*', 'MarkerEdgeColor', '#A2142F'); hold on;
    plot(cost_MOEO(:, 1), cost_MOEO(:, 2), 'd', 'MarkerEdgeColor', '#0072BD'); hold on;
    plot(cost_MOGWEO(:, 1), cost_MOGWEO(:, 2), 'p', 'MarkerEdgeColor', '#EDB120');
    legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel('{\it{F_{loss}}} (kW)', 'FontName', 'Times New Roman', 'FontSize', 16);
    ylabel('{\it{F_{{\Delta}V}}} (%)', 'FontName', 'Times New Roman', 'FontSize', 16);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    grid on;
    axis tight;
elseif M == 3
    figure;
    ax = gca;
    plot3(cost_MOGWO(:, 1), cost_MOGWO(:, 2), cost_MOGWO(:, 3), '*', ...
        'MarkerEdgeColor', '#A2142F'); hold on;
    plot3(cost_MOEO(:, 1), cost_MOEO(:, 2), cost_MOEO(:, 3), 'd', ...
        'MarkerEdgeColor', '#0072BD'); hold on;
    plot3(cost_MOGWEO(:, 1), cost_MOGWEO(:, 2), cost_MOGWEO(:, 3), 'p', ...
        'MarkerEdgeColor', '#EDB120');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', ...
        'FontSize', 16, 'Position', [0.15,0.57,0.27,0.19]);
    xlabel('{\it{F_{loss}}} (kW)', 'FontName', 'Times New Roman', 'FontSize', 16, ...
        'Position', [3.35,2.58,-3.41]);
    ylabel('{\it{F_{{\Delta}V}}} (%)', 'FontName', 'Times New Roman', 'FontSize', 16, ...
        'Position', [2.11,3.04,-3.53]);
    zlabel('{\it{F_{PV}}} (%)', 'FontName', 'Times New Roman', 'FontSize', 16, ...
        'Position', [2.37,3.29,10]);
    grid on;
    set(ax, 'XLim', [2.4, 4], 'XTick', 2.4:0.4:4);
    set(ax, 'YLim', [2.6, 3.2], 'YTick', 2.6:0.15:3.2);
    set(ax, 'ZLim', [0, 20], 'ZTick', 0:5:20);
end
saveas(gca, 'ParetoFront.fig');

toc;

% clear;
% close all;
% clc;
% 
% tic;
% M = 3;
% load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
% Archive_MOGWEO = Archive;
% cost_MOGWEO = reshape([Archive_MOGWEO.Cost], M, [])';
% load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
% Archive_MOGWO = Archive;
% cost_MOGWO = reshape([Archive_MOGWO.Cost], M, [])';
% load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
% Archive_MOEO = Archive;
% cost_MOEO = reshape([Archive_MOEO.Cost], M, [])';
% 
% clear Archive;
% 
% if M == 2
%     figure;
%     ax = gca;
%     ax.FontName = 'Times New Roman';
%     plot(cost_MOGWO(:, 1), cost_MOGWO(:, 2), '*', 'MarkerEdgeColor', '#A2142F'); hold on;
%     plot(cost_MOEO(:, 1), cost_MOEO(:, 2), 'd', 'MarkerEdgeColor', '#0072BD'); hold on;
%     plot(cost_MOGWEO(:, 1), cost_MOGWEO(:, 2), 'p', 'MarkerEdgeColor', '#EDB120');
%     legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', 'FontSize', 10);
%     xlabel('Average system losses (kW)', 'FontName', 'Times New Roman', ...
%         'FontSize', 11);
%     ylabel('Average voltage deviation (%)', 'FontName', 'Times New Roman', ...
%         'FontSize', 11);
%     ax.FontName = 'Times New Roman';
%     ax.FontSize = 11;
%     grid on;
%     axis tight;
%     title('Pareto Optimal Front', 'FontName', 'Times New Roman', 'FontSize', 12);
% elseif M == 3
%     figure;
%     ax = gca;
%     plot3(cost_MOGWO(:, 1), cost_MOGWO(:, 2), cost_MOGWO(:, 3), '*', ...
%         'MarkerEdgeColor', '#A2142F'); hold on;
%     plot3(cost_MOEO(:, 1), cost_MOEO(:, 2), cost_MOEO(:, 3), 'd', ...
%         'MarkerEdgeColor', '#0072BD'); hold on;
%     plot3(cost_MOGWEO(:, 1), cost_MOGWEO(:, 2), cost_MOGWEO(:, 3), 'p', ...
%         'MarkerEdgeColor', '#EDB120');
%     ax.FontName = 'Times New Roman';
%     ax.FontSize = 11;
%     legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', ...
%         'FontSize', 10, 'Position', [0.17,0.63,0.20,0.13]);
%     xlabel('Average system losses (kW)', 'FontName', 'Times New Roman', ...
%         'FontSize', 11);
%     ylabel('Average voltage deviation (%)', 'FontName', 'Times New Roman', ...
%         'FontSize', 11);
%     zlabel('Average curtailment rate (%)', 'FontName', 'Times New Roman', ...
%         'FontSize', 11);
%     grid on;
%     axis tight;
% %     title('Pareto Optimal Front', 'FontName', 'Times New Roman', 'FontSize', 12);
% end
% 
% toc;