clear;
close all;
clc;

tic;
M = 3;
load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWO = Archive;
cost_MOGWO = reshape([Archive_MOGWO.Cost], M, [])';
load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOEO = Archive;
cost_MOEO = reshape([Archive_MOEO.Cost], M, [])';
load('MOGSK_result_of_3objs_under_1_VB.mat', 'pop');
Archive_MOGSK = pop;
cost_MOGSK = reshape([Archive_MOGSK.Cost], M, [])';
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWEO = Archive;
cost_MOGWEO = reshape([Archive_MOGWEO.Cost], M, [])';

min_value = floor(min([min(cost_MOGWO); min(cost_MOEO); ...
    min(cost_MOGSK); min(cost_MOGWEO)], [], 1)*10)/10;
max_value = ceil(max([max(cost_MOGWO); max(cost_MOEO); ...
    min(cost_MOGSK); max(cost_MOGWEO)], [], 1)*10)/10;
interval = 6;

fname = {'{\it{F_{loss}}} (kW)', '{\it{F_{{\Delta}V}}} (%)', '{\it{F_{PV}}} (%)'};

clear Archive;

for count = 1:M
    figure;
    ax = gca;
    h = boxplot([cost_MOGWO(:, count), cost_MOEO(:, count), cost_MOGSK(:, count)...
        cost_MOGWEO(:, count)]);
    set(h, 'LineWidth', 1.5);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    grid on;
    xlabel('Algorithm');
    ylabel(fname{count});
    set(ax, 'XTickLabel', {'MOGWO', 'MOEO', 'MOGSK', 'MOGWEO'}, 'FontSize', 14);
    set(ax, 'YLim', [min_value(count), max_value(count)], ...
        'YTick', linspace(min_value(count), max_value(count), interval));
    saveas(ax, ['Objective Range of F', num2str(count), '.fig']);
end

toc;