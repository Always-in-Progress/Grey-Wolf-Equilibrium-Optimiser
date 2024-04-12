function [Weight, E, Plan] = EWM(Data)
    M = size(Data, 1);
    % Standardisation
    R = (Data - min(Data)) ./ repmat(max(Data)-min(Data), M, 1);
    % Normalisation
    R = R ./ repmat(sum(R), M, 1);
    % Information Entropy Acquisition
    E = -(1/log(M)) * sum(R .* log(R + eps * (R==0)));
    % Attribute Weight Calculation
    Weight = (1 - E) / sum(1-E);
    % Plan Selection
    Plan = sum(R.*repmat(Weight,M,1),2);
end