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
    legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', 'FontSize', 10);
    xlabel('\fontname{宋体}配电网平均损耗\fontname{Times New Roman}(kW)', 'FontSize', 11);
    ylabel('\fontname{宋体}各节点平均电压偏差\fontname{Times New Roman}(%)', 'FontSize', 11);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 11;
    grid on;
    axis tight;
    title('Pareto Optimal Front', 'FontName', 'Times New Roman', 'FontSize', 12);
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
    ax.FontSize = 11;
    legend({'MOGWO', 'MOEO', 'MOGWEO'}, 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Position', [0.17,0.63,0.20,0.13]);
    xlabel('\fontname{宋体}配电网平均损耗\fontname{Times New Roman}(kW)', ...
        'FontSize', 11, 'Position', [3.39,2.60,-1.91]);
    ylabel('\fontname{宋体}各节点平均电压偏差\fontname{Times New Roman}(%)', ...
        'FontSize', 11, 'Position', [2.11,3.07,-1.95]);
    zlabel('\fontname{宋体}分布式光伏平均弃光率\fontname{Times New Roman}(%)',...
        'FontSize', 11, 'Position', [2.33,3.34,10]);
    grid on;
    axis tight;
%     title('Pareto Optimal Front', 'FontName', 'Times New Roman', 'FontSize', 12);
end

toc;