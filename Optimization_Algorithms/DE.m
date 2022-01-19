% The SaNSDE algorithm can be found in:
% Zhenyu Yang, Ke Tang and Xin Yao, "Self-adaptive Differential Evolution with
% Neightborhood Search", in Proceedings of the 2008 IEEE Congress on 
% Evolutionary Computation (CEC2008), Hongkong, China, 2008, pp. 1110-1116.

function [pop, fit, bestsubmem, bestfit, OPTS, FEs]...
    = DE(fname, func_num, pop, gpop, fit, ub, lb, max_iter, dim_index, OPTS)

iter = 0;

[ps, D] = size(pop);

FEs = 0;

while iter < max_iter
iter = iter + 1;

index = zeros(ps, 3);
for i = 1:ps
    index(i, :) = randperm(ps-1, 3);
    index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
end

% 打乱种群
pm1 = pop(index(:, 1),:);             % shuffled population 1
pm2 = pop(index(:, 2),:);             % shuffled population 2
pm3 = pop(index(:, 3),:);             % shuffled population 3

% F and CR
CR = 0.9;
F = 0.5;

%% 变异交叉
aa = rand(ps, D) < CR;
index = find(sum(aa') == 0);
tmpsize = size(index, 2);
for k=1:tmpsize
    bb = ceil(D*rand);
    aa(index(k), bb) = 1;
end
        
mui = aa;
mpo = mui < 0.5;                % inverse mask to mui

ui = pm3 + F .* (pm1 - pm2);
ui = pop .* mpo + ui .* mui;

%% 越界处理
reflect = find(ui > ub);
ui(reflect) = ub - mod((ui(reflect) - ub), (ub - lb));

reflect = find(ui < lb);
ui(reflect) = lb + mod((lb - ui(reflect)), (ub - lb));

gpop(:, dim_index) = ui;

tempfit = feval(fname, gpop, func_num);     % 评估种群适应值
FEs = FEs + ps;

improve = tempfit <= fit;

pop(improve, :) = ui(improve, :);
fit(improve) = tempfit(improve);

end

% get best member and best value
[bestfit, bestid] = min(fit);
bestsubmem = pop(bestid, :);

end
