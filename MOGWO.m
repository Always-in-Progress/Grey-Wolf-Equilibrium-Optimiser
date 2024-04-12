%___________________________________________________________________%
%  Multi-Objective Grey Wolf Optimizer (MOGWO)                      %
%  Source codes demo version 1.0                                    %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%    S. Mirjalili, S. Saremi, S. M. Mirjalili, L. Coelho,           %
%    Multi-objective grey wolf optimizer: A novel algorithm for     %
%    multi-criterion optimization, Expert Systems with Applications,%
%    in press, DOI: http://dx.doi.org/10.1016/j.eswa.2015.10.039    %       %
%                                                                   %
%___________________________________________________________________%

% I acknowledge that this version of MOGWO has been written using
% a large portion of the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

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
M = numel(fopt);

%% GWO Parameters
GreyWolves_num = 200;
MaxIt = 200;  % Maximum Number of Iterations
Archive_size = 200;   % Repository Size

alpha = 0.1;  % Grid Inflation Parameter
nGrid = 10;   % Number of Grids per each Dimension
beta = 4;    % Leader Selection Pressure Parameter
gamma = 2;    % Extra (to be deleted) Repository Member Selection Pressure

% Initialization
GreyWolves = CreateEmptyParticle(GreyWolves_num);

for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        GreyWolves(i).Position(1,j)=unifrnd(min_range(j),max_range(j),1);
    end
    GreyWolves(i).Cost=fobj(GreyWolves(i).Position, fopt, PV);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

costs = GetCosts(GreyWolves);
Archive_costs = GetCosts(Archive);
G = CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

disp(['In iteration ' num2str(0) ': Number of solutions in the archive = ' num2str(numel(Archive))]);

% Results
if drawing_flag
    if M == 2
        hold off;
        plot(costs(1, :), costs(2, :), 'k.');
        hold on;
        plot(Archive_costs(1, :), Archive_costs(2, :), '*');
        legend('Grey wolves','Non-dominated solutions');
        grid on;
        drawnow;
    elseif M == 3
        hold off;
        plot3(costs(1, :), costs(2, :), costs(3, :), 'k.');
        hold on;
        plot3(Archive_costs(1, :), Archive_costs(2, :), Archive_costs(3, :), '*');
        legend('Grey wolves','Non-dominated solutions');
        grid on;
        drawnow;
    end
end

% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    parfor i = 1:GreyWolves_num
        
        % Choose the ALPHA, BETA, and DELTA individuals
        Delta = SelectLeader(Archive,beta);
        Beta = SelectLeader(Archive,beta);
        Alpha = SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive, 1) > 1
            while all(Beta.Position == Delta.Position)
                Beta = SelectLeader(Archive,beta);
            end
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive, 1) > 2
            while all(Alpha.Position == Beta.Position) || all(Alpha.Position == Delta.Position)
                Alpha = SelectLeader(Archive,beta);
            end
        end
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nVar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-A.*abs(D);
        
        % Eq.(3.11) in the paper
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
        % Boundary checking
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,min_range),max_range);
        
        GreyWolves(i).Cost=fobj(GreyWolves(i).Position, fopt, PV);
    end
    
    GreyWolves=DetermineDomination(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);

    % Results
    
    costs = GetCosts(GreyWolves);
    Archive_costs = GetCosts(Archive);
    
    if drawing_flag
        if M == 2
            hold off;
            plot(costs(1, :), costs(2, :), 'k.');
            hold on;
            plot(Archive_costs(1, :), Archive_costs(2, :), '*');
            legend('Grey wolves','Non-dominated solutions');
            grid on;
            drawnow;
        elseif M == 3
            hold off;
            plot3(costs(1, :), costs(2, :), costs(3, :), 'k.');
            hold on;
            plot3(Archive_costs(1, :), Archive_costs(2, :), Archive_costs(3, :), '*');
            legend('Grey wolves','Non-dominated solutions');
            grid on;
            drawnow;
        end
    end
    
end