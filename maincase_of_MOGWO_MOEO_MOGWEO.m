clear;
close all;
clc;

fprintf('The program starts running at '); disp(datetime(clock));
tic;

disp('MOGWO is employed to solve the problem');
MOGWO;
save('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat');

fprintf('\n');
disp('MOEO is employed to solve the problem');
MOEO;
save('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat');

fprintf('\n');
disp('MOGWEO is employed to solve the problem');
MOGWEO;
save('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat');

toc;

fprintf('\n');
fprintf('The program ends running at '); disp(datetime(clock));