function [w, CI, CR] = AHP(A, err)
% EIGVET ����ֵ�����������������������
    if ~exist('err', 'var')
        err = 1e-3;
    end
    % �ж��Ƿ�Ϊ����
    [minput, ninput] = size(A);
    if minput ~= ninput
        disp('���Ƿ���, ����������!');
        return;
    end
    
    % �ж��Ƿ�Ϊ����������
    t = A .* (A.');
    if ~all(A(:) > 0) || ~all(abs(t(:) - 1) <= err)
        disp('��������������, ����������!');
        return;
    end
    
    [W, lamda] = eig(A);
    [~, n] = find(lamda == max(max(lamda)));
    w = W(:, n) ./ sum(W(:, n));

    lamda = max(max(eig(A)));
    CI = (lamda-length(A)) / (length(A) - 1);
    CR = CI / ri(length(A));
    
end