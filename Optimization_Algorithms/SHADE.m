function [pop, fit, subbestmem, bestfit, OPTS, used_FEs, current_pop, current_fit, success_rate] = SHADE(fname, func_num, pop, gpop, fit, ub, lb, max_iter, dim_index, OPTS)


[NP, D] = size(pop);
% 算法初始参数
H = 6;
pmin = 2/NP;
Max_FEs = 3e6;

% rNinit = 18;
rarc = 1;
% Ninit = round(D * rNinit);
Ninit = 100;

if (OPTS.first == 1)
    A_size = round(Ninit * rarc);
    
    MCR = 0.5 * ones(H, 1);
    MF = 0.5 * ones(H, 1);
    A = [];     % archive
    
%     iters = 1;
    k = 1;
    
    FEs = 0;
    
    OPTS.first = 0;
else
    MCR     = OPTS.MCR;
    MF      = OPTS.MF;
%     iters   = OPTS.iters;
    k       = OPTS.k;
    A       = OPTS.A;
    A_size  = OPTS.A_size;
    FEs     = OPTS.FEs;
end

[bestfit, ibest] = min(fit);
subbestmem = pop(ibest, :);
used_FEs = 0;

current_pop = num2cell(pop, 2);
current_fit = num2cell(fit);
success_num = 0;

% terminal_iters = iters + max_iters;
ava_FEs = 100 * max_iter;
iters = 0;

while used_FEs < ava_FEs
    SCR = [];
    SF = [];
    delta_fit = [];
    
    % 计算 CR 和 F
    r = randi(H, 1, NP);
    
    CR_sigma = MCR(r);
    CR = normrnd(CR_sigma, 0.1, NP, 1);
    CR(CR == Inf) = 0;
    CR(CR > 1) = 1;
    CR(CR < 0) = 0;
    
    F_sigma = MF(r);
    F = normrnd(F_sigma, 0.1, NP, 1);
    F(F > 1) = 1;
    
    index = find(F < 0);
    for i = 1 : length(index)
        F(index(i)) = normrnd(F_sigma(i), 0.1);
        while F(index(i)) < 0
            F(index(i)) = normrnd(F_sigma(i), 0.1);
        end
    end
    
    % 随机数小于CR为1，否则为0
    aa = rand(NP, D) < repmat(CR, 1, D);
    index = find(sum(aa') == 0);
    tmpsize = size(index, 2);
    for t = 1 : tmpsize
        bb = ceil(D * rand);
        aa(index(t), bb) = 1;
    end
    mui = aa;
    mpo = mui < 0.5;
    
    % 计算pbest, pm1, pm2
    index = zeros(NP, 2);
    index(:, 1) = randi(NP - 1, NP, 1);
    index(:, 2) = randi(NP + size(A, 1) - 2, NP, 1);
    for i = 1 : NP
        index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
        if index(i, 2) >= index(i, 1)
            index(i, 2) = index(i, 2) + 1;
        end
    end
    
    pm1 = pop(index(:, 1), :);             % shuffled population 1
    
    union_pop = [pop; A];
    pm2 = union_pop(index(:, 2), :);       % shuffled population 2
    
    [~, sidx] = sort(fit);
    p = pmin + rand(NP, 1) * (0.2 - pmin);
    pbest = zeros(NP, D);
    for i = 1 : NP
        pbest_index = sidx(randi(ceil(p(i) * NP)));
        pbest(i, :) = pop(pbest_index, :); 
    end

    % 交叉变异 current-to-pbest/1/bin
    ui = pop + repmat(F, 1, D) .* (pbest - pop) + repmat(F, 1, D) .* (pm1 - pm2);
    ui = pop .* mpo + ui .* mui;
    
    % 越界处理
    repair = find(ui > ub);
    ui(repair) = (ub + pop(repair))/2;
    repair = find(ui < lb);
    ui(repair) = (lb + pop(repair))/2;
    
    % 评估
    gpop(:, dim_index) = ui;
    tempfit = feval(fname, gpop, func_num);
    used_FEs = used_FEs + NP;
    
    % 选择
    is_improved = tempfit < fit;
    A = [A; pop(is_improved, :)];
    SCR = [SCR; CR(is_improved)];
    SF = [SF; F(is_improved)];
    delta_fit = [delta_fit, abs(fit(is_improved) - tempfit(is_improved))];
    
    is_updated = tempfit <= fit;
    pop(is_updated, :) = ui(is_updated, :);
    fit(is_updated) = tempfit(is_updated);
    
    for i = 1 : NP
        current_pop{i} = [current_pop{i}; ui(i, :)];
        current_fit{i} = [current_fit{i}, tempfit(i)];
    end
    
    % 调整记忆集A的大小
    if size(A, 1) > A_size
        delete_num = size(A, 1) - A_size;
        A(randperm(size(A, 1), delete_num), :) = [];
    end
    
    % 更新最优解和最优个体
    [temp_bestfit, ibest] = min(fit);
    temp_subbestmem = pop(ibest, :);
    
    if temp_bestfit < bestfit
        bestfit = temp_bestfit;
        subbestmem = temp_subbestmem;
        success_num = success_num + 1;
    end
    
    % fprintf('Gen: %d, Func Val: %f\n', iters, bestfit - temp_bestfit);
    
    % 更新MCR和MF
    [MCR, MF, k] = UpdateMemory(MCR, MF, SCR, SF, delta_fit, k);
    
%     % 种群规模线性下降
%     N = round((Nmin - Ninit)/Max_FEs * (FEs+used_FEs) + Ninit);
%     if N < NP
%         [~, sidx] = sort(fit);
%         pop = pop(sidx(1:N), :);
%         fit = fit(sidx(1:N));
%         
%         A_size = round(N * rarc);
%     end
    
%     NP = size(pop, 1);
%     gpop = gpop(1:NP, :);
    
    iters = iters + 1;
end

gpop(:, dim_index) = pop; 
success_rate = success_num ./ iters;


OPTS.MCR     = MCR;
OPTS.MF      = MF;
OPTS.iters   = iters;
OPTS.k       = k;
OPTS.A       = A;
OPTS.A_size  = A_size;

end

function [MCR, MF, k] = UpdateMemory(MCR, MF, SCR, SF, delta_fit, k)
    w = delta_fit ./ sum(delta_fit);
    if ~isempty(SCR) && ~isempty(SF)
        if MCR(k) == Inf || max(SCR) == 0
            MCR(k) = Inf;
        else
            MCR(k) = (w * SCR .^ 2) / (w * SCR);
        end
        MF(k) = (w * SF .^ 2) / (w * SF);
        k = k + 1;
        if k > length(MCR)
            k = 1;
        end
    end
end
