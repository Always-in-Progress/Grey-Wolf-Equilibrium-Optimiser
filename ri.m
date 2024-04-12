function Rri = ri( num )
% Rri = ri( num )
% 求num阶正互反一致性判断矩阵的随机一致性指标值

    % 初始化
    times = 1000;
    Lamda = zeros(1,times);
    ci = zeros(1,times);
    for num1 = 1:length(num)
        for count = 1:times
            A = ones(num(num1));
            for k = 1:num(num1)
                for l = k:num(num1)
                    if k ~= l
                        flag = round(rand(1));
                        intr = round(8*rand(1))+1;
                        A(k,l) = flag*intr+(1-flag)*1/intr;
                        A(l,k) = 1/A(k,l);
                    end
                end
            end
            [alpha,lamda] = eig(A);
            Lamda(count) = max(max(lamda));
            if num(num1) > 2
                ci(count) = (Lamda(count)-num(num1))/(num(num1)-1);
            else
                if Lamda(count) == num(num1)
                    ci(count) = 0;
                end
            end
        end
        Rri(num1) = sum(ci)/times;
    end
end