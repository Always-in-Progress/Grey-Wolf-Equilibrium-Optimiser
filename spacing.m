clear;
close all;
clc;

tic;

%% SP of MOGWO
load('MOGWO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'M');

Archive_size = numel(Archive);
cost = zeros(Archive_size, M);

for count = 1:Archive_size
    cost(count, :) = Archive(count).Cost;
end

% Normalization
uni_cost = (cost - min(cost)) ./ (max(cost) - min(cost));

distance = zeros(Archive_size, Archive_size);

for row = 1:Archive_size
    for col = 1:Archive_size
        if row ~= col
            distance(row, col) = sqrt(sum((uni_cost(row, :) - uni_cost(col, :)).^2));
        else
            distance(row, col) = inf;
        end
    end
end

d = min(distance, [], 2);

SP_MOGWO = std(d);

%% SP of MOEO
load('MOEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'M');

Archive_size = numel(Archive);
cost = zeros(Archive_size, M);

for count = 1:Archive_size
    cost(count, :) = Archive(count).Cost;
end

% Normalization
uni_cost = (cost - min(cost)) ./ (max(cost) - min(cost));

distance = zeros(Archive_size, Archive_size);

for row = 1:Archive_size
    for col = 1:Archive_size
        if row ~= col
            distance(row, col) = sqrt(sum((uni_cost(row, :) - uni_cost(col, :)).^2));
        else
            distance(row, col) = inf;
        end
    end
end

d = min(distance, [], 2);

SP_MOEO = std(d);

%% SP of MOGSK
load('MOGSK_result_of_3objs_under_1_VB.mat', 'pop', 'nObj');

Archive_size = numel(pop);
cost = zeros(Archive_size, nObj);

for count = 1:Archive_size
    cost(count, :) = pop(count).Cost;
end

% Normalization
uni_cost = (cost - min(cost)) ./ (max(cost) - min(cost));

distance = zeros(Archive_size, Archive_size);

for row = 1:Archive_size
    for col = 1:Archive_size
        if row ~= col
            distance(row, col) = sqrt(sum((uni_cost(row, :) - uni_cost(col, :)).^2));
        else
            distance(row, col) = inf;
        end
    end
end

d = min(distance, [], 2);

SP_MOGSK = std(d);

%% SP of MOGWEO
load('MOGWEO_result_of_3objs_under_1_VB_in_10_grids.mat', 'Archive', 'M');

Archive_size = numel(Archive);
cost = zeros(Archive_size, M);

for count = 1:Archive_size
    cost(count, :) = Archive(count).Cost;
end

% Normalization
uni_cost = (cost - min(cost)) ./ (max(cost) - min(cost));

distance = zeros(Archive_size, Archive_size);

for row = 1:Archive_size
    for col = 1:Archive_size
        if row ~= col
            distance(row, col) = sqrt(sum((uni_cost(row, :) - uni_cost(col, :)).^2));
        else
            distance(row, col) = inf;
        end
    end
end

d = min(distance, [], 2);

SP_MOGWEO = std(d);

save('spacing.mat', 'SP_MOGWO', 'SP_MOEO', 'SP_MOGSK', 'SP_MOGWEO');

toc;