function [w, CI, CR] = AHP(A, err)
% EIGVET 特征值法求正互反矩阵的排序向量
    if ~exist('err', 'var')
        err = 1e-3;
    end
    % 判断是否为方阵
    [minput, ninput] = size(A);
    if minput ~= ninput
        disp('不是方阵, 请重新输入!');
        return;
    end
    
    % 判断是否为正互反矩阵
    t = A .* (A.');
    if ~all(A(:) > 0) || ~all(abs(t(:) - 1) <= err)
        disp('不是正互反矩阵, 请重新输入!');
        return;
    end
    
    [W, lamda] = eig(A);
    [~, n] = find(lamda == max(max(lamda)));
    w = W(:, n) ./ sum(W(:, n));

    lamda = max(max(eig(A)));
    CI = (lamda-length(A)) / (length(A) - 1);
    CR = CI / ri(length(A));
    
end