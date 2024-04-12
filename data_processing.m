clear;
close all;
clc;

tic;
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'fopt', ...
    'PV', 'nVar', 'M');
% load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'fopt', ...
%     'PV', 'nVar', 'M');
% load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'fopt', ...
%     'PV', 'nVar', 'M');
% load('MOGWEO_result_of_3objs_under_1_VB_in_30_grids.mat', 'Archive', 'fopt', ...
%     'PV', 'nVar', 'M');
% load('MOGWO_result_of_3objs_under_1_VB_in_30_grids.mat', 'Archive', 'fopt', ...
%     'PV', 'nVar', 'M');
% load('MOEO_result_of_3objs_under_1_VB_in_30_grids.mat', 'Archive', 'fopt', ...
%     'PV', 'nVar', 'M');

draw_flag = true;

Archive_size = numel(Archive);
sol = reshape(struct, Archive_size, 1);
sort_cost = zeros(Archive_size, M);

for n = 1:Archive_size
    x = Archive(n).Position;
    [f, ploss, vol_pu, vol_deg, P, Q] = fitness(x, fopt, PV);
    sol(n).allocation = decouple(x);
    sol(n).fcost = f;
    sol(n).ploss = ploss;
    sol(n).vol_pu = vol_pu;
    sol(n).vol_deg = vol_deg;
    sol(n).P = P;
    sol(n).Q = Q;
    sort_cost(n, :) = sol(n).fcost;
end

% 归一化
uni_sort_cost = (sort_cost - min(sort_cost)) ./ (max(sort_cost) - min(sort_cost));

[~, index1] = min(sort_cost(:, 1));
[~, index2] = min(sort_cost(:, 2));
[~, index3] = min(sort_cost(:, 3));
[~, index4] = min(sqrt(sum(uni_sort_cost.^2, 2)));

[f_0, ploss_0, vol_pu_0, vol_deg_0, P_0, Q_0] = fitness(zeros(1, nVar), fopt, PV);

cost_of_dispatch = [f_0; sort_cost([index1, index2, index3, index4], :);];

if draw_flag    
    if M == 2
        figure;
        plot(sort_cost(:, 1), sort_cost(:, 2), '*');
        xlabel('配电网平均损耗(kW)', 'FontName', '宋体');
        ylabel('各节点平均电压偏差(%)', 'FontName', '宋体');
        grid on;
        axis tight;
        title('Pareto Optimal Front');
    elseif M == 3
        figure;
        plot3(sort_cost(:, 1), sort_cost(:, 2), sort_cost(:, 3), '*');
        xlabel('配电网平均损耗(kW)', 'FontName', '宋体');
        ylabel('各节点平均电压偏差(%)', 'FontName', '宋体');
        zlabel('分布式光伏平均弃光率(%)', 'FontName', '宋体');
        grid on;
        axis tight;
        title('Pareto Optimal Surface');
    end
    
    figure;
    bar3(0:23, ploss_0*1000);
    axis tight;
    xlabel('支路编号', 'FontName', 'Times New Roman + SimSung');
    ylabel('时段(h)', 'FontName', 'Times New Roman + SimSung');
    zlabel('支路损耗(kW)', 'FontName', 'Times New Roman + SimSung');
    title('不进行补偿优化时不同时段各支路的损耗情况', 'FontName', 'Times New Roman + SimSung');
    figure;
    bar3(0:23, sol(index1).ploss*1000);
    axis tight;
    xlabel('支路编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('支路损耗(kW)', 'FontName', '宋体');
    title('利用单因素补偿优化方案1优化时不同时段各支路的损耗情况', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index2).ploss*1000);
    axis tight;
    xlabel('支路编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('支路损耗(kW)', 'FontName', '宋体');
    title('利用单因素补偿优化方案2优化时不同时段各支路的损耗情况', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index3).ploss*1000);
    axis tight;
    xlabel('支路编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('支路损耗(kW)', 'FontName', '宋体');
    title('利用单因素补偿优化方案3优化时不同时段各支路的损耗情况', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index4).ploss*1000);
    axis tight;
    xlabel('支路编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('支路损耗(kW)', 'FontName', '宋体');
    title('利用多因素折中补偿优化方案优化时不同时段各支路的损耗情况', 'FontName', '宋体');
    
    figure;
    plot(0:23, sum(ploss_0, 2)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sum(sol(index1).ploss, 2)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sum(sol(index2).ploss, 2)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sum(sol(index3).ploss, 2)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sum(sol(index4).ploss, 2)*1000, 'o-', 'LineWidth', 1.2);
    grid on;
    axis tight;
    legend({'不进行补偿优化', '利用单因素补偿优化方案1进行补偿优化', '利用单因素补偿优化方案2进行补偿优化', ...
        '利用单因素补偿优化方案3进行补偿优化', '利用多因素折中补偿优化方案进行补偿优化'}, 'Location', ...
        'NorthWest', 'FontName', '宋体');
    xlabel('时段(h)', 'FontName', '宋体');
    ylabel('总系统损耗(kW)', 'FontName', '宋体');
    
    figure;
    bar3(0:23, vol_pu_0);
    zlim([0.95 1.07]);
    xlabel('节点编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('电压标幺值(p.u.)', 'FontName', '宋体');
    title('不进行补偿优化时不同时段各母线的电压分布', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index1).vol_pu);
    zlim([0.95 1.07]);
    xlabel('节点编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('电压标幺值(p.u.)', 'FontName', '宋体');
    title('利用单因素补偿优化方案1进行补偿优化时不同时段各母线的电压分布', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index2).vol_pu);
    zlim([0.95 1.07]);
    xlabel('节点编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('电压标幺值(p.u.)', 'FontName', '宋体');
    title('利用单因素补偿优化方案2进行补偿优化时不同时段各母线的电压分布', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index3).vol_pu);
    zlim([0.95 1.07]);
    xlabel('节点编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('电压标幺值(p.u.)', 'FontName', '宋体');
    title('利用单因素补偿优化方案3进行补偿优化时不同时段各母线的电压分布', 'FontName', '宋体');
    figure;
    bar3(0:23, sol(index4).vol_pu);
    zlim([0.95 1.07]);
    xlabel('节点编号', 'FontName', '宋体');
    ylabel('时段(h)', 'FontName', '宋体');
    zlabel('电压标幺值(p.u.)', 'FontName', '宋体');
    title('利用多因素折中补偿优化方案进行补偿优化时不同时段各母线的电压分布', 'FontName', '宋体');
    
    figure;
    plot(0:23, mean(abs(vol_pu_0 - 1), 2)*100, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, mean(abs(sol(index1).vol_pu - 1), 2)*100, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, mean(abs(sol(index2).vol_pu - 1), 2)*100, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, mean(abs(sol(index3).vol_pu - 1), 2)*100, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, mean(abs(sol(index4).vol_pu - 1), 2)*100, 'o-', 'LineWidth', 1.2);
    grid on;
    axis tight;
    legend({'不进行补偿优化', '利用单因素补偿优化方案1进行补偿优化', '利用单因素补偿优化方案2进行补偿优化', ...
        '利用多因素折中补偿优化方案进行补偿优化'}, 'Location', 'NorthWest', ...
        'FontName', '宋体');
    xlabel('时段(h)', 'FontName', '宋体');
    ylabel('平均电压偏差(%)', 'FontName', '宋体');
    
    figure;
    plot(0:23, PV*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sol(index1).allocation(:, 13)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sol(index2).allocation(:, 13)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sol(index3).allocation(:, 13)*1000, 'o-', 'LineWidth', 1.2); hold on;
    plot(0:23, sol(index4).allocation(:, 13)*1000, 'o-', 'LineWidth', 1.2);
    grid on;
    axis tight;
    legend({'分布式光伏最大出力', '单因素补偿优化方案1对应的分布式光伏出力', '单因素补偿优化方案2对应的分布式光伏出力', ...
        '单因素补偿优化方案3对应的分布式光伏出力', '多因素折中补偿优化方案的分布式光伏出力'}, 'Location', 'NorthWest', ...
        'FontName', '宋体');
    xlabel('时段(h)', 'FontName', '宋体');
    ylabel('分布式光伏出力(kW)', 'FontName', '宋体');
end
toc;