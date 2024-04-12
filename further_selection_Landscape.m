clear;
close all;
clc;

nMethod = 10;
nObj = 3;
nAlgo = 4;

%% Preliminaries
% Retrieve cost data
% MOGWO
load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWO = Archive;
cost_MOGWO = reshape([Archive_MOGWO.Cost], nObj, [])';
% MOEO
load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOEO = Archive;
cost_MOEO = reshape([Archive_MOEO.Cost], nObj, [])';
% MOGSK
load('MOGSK_result_of_3objs_under_1_VB', 'pop');
Archive_MOGSK = pop;
cost_MOGSK = reshape([Archive_MOGSK.Cost], nObj, [])';
% MOGWEO
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive');
Archive_MOGWEO = Archive;
cost_MOGWEO = reshape([Archive_MOGWEO.Cost], nObj, [])';

%% Normalisation
% Pre-calculation
max_cost = max([cost_MOGWO; cost_MOEO; cost_MOGSK; cost_MOGWEO], [], 1);
min_cost = min([cost_MOGWO; cost_MOEO; cost_MOGSK; cost_MOGWEO], [], 1);
% Transform the cost into a normalised value
normalized_cost_MOGWO = (cost_MOGWO - min_cost) ./ (max_cost - min_cost);
normalized_cost_MOEO = (cost_MOEO - min_cost) ./ (max_cost - min_cost);
normalized_cost_MOGSK = (cost_MOGSK - min_cost) ./ (max_cost - min_cost);
normalized_cost_MOGWEO = (cost_MOGWEO - min_cost) ./ (max_cost - min_cost);

%% Weight Acquisition
Weight = zeros(nMethod, nObj);
% Method 1-3: Only one objective function is considered
% Method 1: Only F_loss is considered
Weight(1, :) = [1, 0, 0];
% Method 2: Only F_ΔV is considered
Weight(2, :) = [0, 1, 0];
% Method 3: Only F_PV is considered
Weight(3, :) = [0, 0, 1];
% Method 4-9: The attribute weight of each objective function is obtained by AHP
Weight(4, :) = AHP(ones(nObj, nObj));
Weight(5, :) = AHP([1, 1/2, 1/3; 2, 1, 1; 3, 1, 1]);
Weight(6, :) = AHP([1, 1, 5; 1, 1, 5; 1/5, 1/5, 1]);
Weight(7, :) = AHP([1, 4.327, 0.481; 0.231, 1, 0.1111; 2.08, 9, 1]);
Weight(8, :) = AHP([1, 1/3, 2; 3, 1, 6; 1/2, 1/6, 1]);
Weight(9, :) = AHP([1, 8, 2; 1/8, 1, 1/5; 1/2, 5, 1]);
% Method 10: The attribute weight of each objective function is acquired by EWM
Weight(10, :) = EWM(cost_MOGWEO);

%% Plan Selection
% Comprehensive cost acquisition
comprehensive_cost_MOGWO = normalized_cost_MOGWO * Weight';
comprehensive_cost_MOEO = normalized_cost_MOEO * Weight';
comprehensive_cost_MOGSK = normalized_cost_MOGSK * Weight';
comprehensive_cost_MOGWEO = normalized_cost_MOGWEO * Weight';
% Obtain the Minimum Comprehensive Cost
% MOGWO
[min_comprehensive_cost_MOGWO, flag_MOGWO] = min(comprehensive_cost_MOGWO, [], 1);
min_comprehensive_cost_MOGWO = min_comprehensive_cost_MOGWO';
flag_MOGWO = flag_MOGWO';
% MOEO
[min_comprehensive_cost_MOEO, flag_MOEO] = min(comprehensive_cost_MOEO, [], 1);
min_comprehensive_cost_MOEO = min_comprehensive_cost_MOEO';
flag_MOEO = flag_MOEO';
% MOGSK
[min_comprehensive_cost_MOGSK, flag_MOGSK] = min(comprehensive_cost_MOGSK, [], 1);
min_comprehensive_cost_MOGSK = min_comprehensive_cost_MOGSK';
flag_MOGSK = flag_MOGSK';
% MOGWEO
[min_comprehensive_cost_MOGWEO, flag_MOGWEO] = min(comprehensive_cost_MOGWEO, [], 1);
min_comprehensive_cost_MOGWEO = min_comprehensive_cost_MOGWEO';
flag_MOGWEO = flag_MOGWEO';
% Gather overall information
info = [Weight, cost_MOGWO(flag_MOGWO, :), min_comprehensive_cost_MOGWO, ...
    cost_MOEO(flag_MOEO, :), min_comprehensive_cost_MOEO, ...
    cost_MOGSK(flag_MOGSK, :), min_comprehensive_cost_MOGSK, ...
    cost_MOGWEO(flag_MOGWEO, :), min_comprehensive_cost_MOGWEO];
info_format = round(info * 1e4) / 1e4;
% Improvement of F^*
improF = (info(:, 7:4:15) - repmat(info(:, end), 1, nAlgo - 1)) ./ ...
    info(:, 7:4:15) * 100;
improF_format = round(improF * 1e4) / 1e4;