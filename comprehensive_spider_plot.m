clear;
close all;
clc;

M = 3;

%% Preliminaries
% Retrieve cost data
% MOGWO
min_MOGWO = [2.6002, 2.7254, 2.7439];
load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWO = Archive;
cost_MOGWO = reshape([Archive_MOGWO.Cost], M, [])';
range_MOGWO = max(cost_MOGWO, [], 1) - min(cost_MOGWO, [], 1);
% MOEO
min_MOEO = [2.3939, 2.6773, 0.0130];
load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOEO = Archive;
cost_MOEO = reshape([Archive_MOEO.Cost], M, [])';
range_MOEO = max(cost_MOEO, [], 1) - min(cost_MOEO, [], 1);
% MOGSK
min_MOGSK = [2.5453, 2.9246, 5.1766];
load('MOGSK_result_of_3objs_under_1_VB', 'pop');
Archive_MOGSK = pop;
cost_MOGSK = reshape([Archive_MOGSK.Cost], M, [])';
range_MOGSK = max(cost_MOGSK, [], 1) - min(cost_MOGSK, [], 1);
% MOGWEO
min_MOGWEO = [2.3782, 2.6861, 0];
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWEO = Archive;
cost_MOGWEO = reshape([Archive_MOGWEO.Cost], M, [])';
range_MOGWEO = max(cost_MOGWEO, [], 1) - min(cost_MOGWEO, [], 1);

load('spacing.mat');
data = [min_MOGWO, range_MOGWO, SP_MOGWO
    min_MOEO, range_MOEO, SP_MOEO
    min_MOGSK, range_MOGSK, SP_MOGSK
    min_MOGWEO, range_MOGWEO, SP_MOGWEO];

%% Spider chart in English Version
figure;
spider_plot(data, 'FillOption', 'on', ...
    'AxesLabels', {'min {\it{F_{loss}}} (kW)', 'min {\it{F_{\DeltaV}}} (%)', ...
    'min {\it{F_{PV}}} (%)', 'Range of {\it{F_{loss}}} (kW)', ...
    'Range of {\it{F_{\DeltaV}}} (%)', 'Range of {\it{F_{PV}}} (%)', ...
    '{\it{SP}}'}, ...
    'FillTransparency', 0.2, ...
    'AxesDirection', {'reverse', 'reverse', 'reverse', 'normal', 'normal', ...
    'normal', 'reverse'}, ...
    'AxesLabelsEdge', 'none', ...
    'Color', [0, 114, 189; 217, 83, 25; 126, 47, 142; 119, 172, 48]/255,...
    'AxesScaling', 'linear', ...
    'AxesPrecision', [2, 2, 2, 2, 2, 2, 4], ...
    'AxesFont', 'Times New Roman', ...
    'LabelFont', 'Times New Roman', ...
    'AxesFontSize', 14, ...
    'LabelFontSize', 16);
legend({'GWO', 'EO', 'GSK', 'GWEO'}, 'FontName', 'Times New Roman',...
    'FontSize', 16, 'Position', [0.74, 0.8, 0.19, 0.19]);
set(gcf, 'Position', [350, 160, 762, 552]);