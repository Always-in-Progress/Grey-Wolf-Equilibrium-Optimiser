%% Optimal Power Flow Utilising MOGSK Based on NSGA-?
clear;
close all;
clc;

tic;
mpc = case118zh;

% 粒子群算法(PSO)、樽海鞘群算法(SSA)、灰色狼群算法(GWO)结果的综合比较
N = 48; % 补偿装置总数
nVar = 24*(N+1); %维度（决策变量的个数）
% 无功补偿设备
QL = sort(mpc.bus(:, 4), 'descend')';
QC = QL(1:N);
PV = [0,0,0,0,0,0,0.0670217995175985,0.204959666510439,0.441740303387099,0.674386614692288,0.857661076753659,0.975895167221085,0.815624321265675,1,0.908416115105721,0.633496801137035,0.387087186639079,0.317973456486700,0.0825036103903032,0,0,0,0,0]';
PV = PV * 0.5;

drawing_flag = true;

% 下边界
Lb = zeros(24, N);
Lb(:, N+1) = 0.8*PV;
% 上边界
Ub = repmat(QC, 24, 1);
Ub(:, N+1) = PV;
% 改变维度布局
min_range = reshape(Lb', 1, nVar);
max_range = reshape(Ub', 1, nVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective Functions
fopt = {'mean(mean(ploss)) * 1000', 'mean(mean(abs(vol_pu - 1))) * 100', ...
    '(1 - mean(P_PV(PV > 0)./PV(PV > 0))) * 100'};
fobj = @fitness;
% Number of Objective Functions
nObj = numel(fopt);

%% MOGSK Parameters
MaxIt = 200;  % Maximum Number of Iterations
nPop = 200;  % Population Size
max_nfes = 10000 * nVar;

G_Max = fix(max_nfes/nPop);
nfes = 0;

KF = 0.5;% Knowledge Factor
KR = 0.9;%Knowledge Ratio
K = 10*ones(nPop, 1);%Knowledge Rate
lu = [min_range; max_range];

%% Initialization
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    for j = 1:nVar
        pop(i, 1).Position(1, j) = unifrnd(min_range(j),max_range(j),1);
    end
    pop(i, 1).Cost = fobj(pop(i).Position, fopt, PV);
end

% Non-Dominated Sorting
[pop, F] = NonDominatedSorting(pop);

% Calculate Crowding Distance
pop = CalcCrowdingDistance(pop, F);

% Sort Population
[pop, F] = SortPopulation(pop);

% Store F1
F1 = pop(F{1});

disp('Staring MOGSK ...');
% Show Iteration Information
disp(['Iteration ' num2str(0) ': Number of F1 Members = ' num2str(numel(F1))]);

popCosts = reshape([pop.Cost], nObj, nPop);
popRanks = [pop.Rank];

% Results
if drawing_flag
    if nObj == 2
        hold off;
        plot(popCosts(1, :), popCosts(2, :), 'k.');
        hold on;
        plot(popCosts(1, popRanks == 1), popCosts(2, popRanks == 1), '*');
        legend('Individuals', 'Non-dominated solutions');
        grid on;
        drawnow;
    elseif nObj == 3
        hold off;
        plot3(popCosts(1, :), popCosts(2, :), popCosts(3, :), 'k.');
        hold on;
        plot3(popCosts(1, popRanks == 1), popCosts(2, popRanks == 1), ...
            popCosts(3, popRanks == 1), '*');
        legend('Individuals', 'Non-dominated solutions');
        grid on;
        drawnow;
    end
end

% Iteration
for iter = 1:MaxIt
    D_Gained_Shared_Junior = ceil((nVar)*(1-iter/G_Max).^K);
    D_Gained_Shared_Senior = nVar-D_Gained_Shared_Junior;
    %     pop = pop; % the old population becomes the current population
    pos = reshape([pop.Position]', nVar, nPop)';

    indBest = 1:nPop;
    [Rg1, Rg2, Rg3] = Gained_Shared_Junior_R1R2R3(indBest);

    [R1, R2, R3] = Gained_Shared_Senior_R1R2R3(indBest);
    R01 = 1:nPop;
    Gained_Shared_Junior=zeros(nPop, nVar);

    ind1 = R01 > Rg3;
    if(sum(ind1) > 0)
        Gained_Shared_Junior(ind1,:)= pos(ind1,:) + KF*ones(sum(ind1), nVar) .* (pos(Rg1(ind1),:) - pos(Rg2(ind1),:)+pos(Rg3(ind1), :)-pos(ind1,:)) ;
    end
    ind1 = ~ind1;
    if(sum(ind1) > 0)
        Gained_Shared_Junior(ind1,:) = pos(ind1,:) + KF*ones(sum(ind1), nVar) .* (pos(Rg1(ind1),:) - pos(Rg2(ind1),:)+pos(ind1,:)-pos(Rg3(ind1), :)) ;
    end
    R0 = 1:nPop;
    Gained_Shared_Senior = zeros(nPop, nVar);
    ind = R0 > R2;
    if(sum(ind) > 0)
        Gained_Shared_Senior(ind,:) = pos(ind,:) + KF*ones(sum(ind), nVar) .* (pos(R1(ind),:) - pos(ind,:) + pos(R2(ind),:) - pos(R3(ind), :)) ;
    end
    ind = ~ind;
    if(sum(ind)>0)
        Gained_Shared_Senior(ind,:) = pos(ind,:) + KF*ones(sum(ind), nVar) .* (pos(R1(ind),:) - pos(R2(ind),:) + pos(ind,:) - pos(R3(ind), :)) ;
    end
    Gained_Shared_Junior = boundConstraint(Gained_Shared_Junior, pos, lu);
    Gained_Shared_Senior = boundConstraint(Gained_Shared_Senior, pos, lu);

    D_Gained_Shared_Junior_mask = rand(nPop, nVar) <= (D_Gained_Shared_Junior(:, ones(1, nVar))./nVar);
    D_Gained_Shared_Senior_mask = ~D_Gained_Shared_Junior_mask;

    D_Gained_Shared_Junior_rand_mask = rand(nPop, nVar) <= KR*ones(nPop, nVar);
    D_Gained_Shared_Junior_mask = and(D_Gained_Shared_Junior_mask, D_Gained_Shared_Junior_rand_mask);

    D_Gained_Shared_Senior_rand_mask = rand(nPop, nVar)<=KR*ones(nPop, nVar);
    D_Gained_Shared_Senior_mask = and(D_Gained_Shared_Senior_mask, D_Gained_Shared_Senior_rand_mask);
    ui = pos;

    ui(D_Gained_Shared_Junior_mask) = Gained_Shared_Junior(D_Gained_Shared_Junior_mask);
    ui(D_Gained_Shared_Senior_mask) = Gained_Shared_Senior(D_Gained_Shared_Senior_mask);

    popnew = repmat(empty_individual, nPop, 1);
    parfor k = 1:nPop
        popnew(k, 1).Position = ui(k, :);
        popnew(k, 1).Cost = fobj(popnew(k, 1).Position, fopt, PV);
    end

    % Merge
    pop = [pop; popnew];

    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    pop = SortPopulation(pop);

    % Truncate
    pop = pop(1:nPop);

    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    [pop, F] = SortPopulation(pop);

    % Store F1
    F1 = pop(F{1});

    % Show Iteration Information
    disp(['Iteration ' num2str(iter) ': Number of F1 Members = ' num2str(numel(F1))]);

    popCosts = reshape([pop.Cost], nObj, nPop);
    popRanks = [pop.Rank];

    % Results
    if drawing_flag
        if nObj == 2
            hold off;
            plot(popCosts(1, :), popCosts(2, :), 'k.');
            hold on;
            plot(popCosts(1, popRanks == 1), popCosts(2, popRanks == 1), '*');
            legend('Individuals', 'Non-dominated solutions');
            grid on;
            drawnow;
        elseif nObj == 3
            hold off;
            plot3(popCosts(1, :), popCosts(2, :), popCosts(3, :), 'k.');
            hold on;
            plot3(popCosts(1, popRanks == 1), popCosts(2, popRanks == 1), ...
                popCosts(3, popRanks == 1), '*');
            legend('Individuals', 'Non-dominated solutions');
            grid on;
            drawnow;
        end
    end
end

toc;