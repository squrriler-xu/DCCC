% The DCCC algorithm can be found in:
% Peilan Xu, Wenjian Luo, Xin Lin, Yatong Chang, Yuhui Shi.
% "Cooperative Coevolution Based on the Difficulty and Contribution for Large-Scale Optimization"
% IEEE Transcations on Evolutionary Computation, Revised and ReSubmitted

function [bestfit] = DCCC_final(benchmark, fname, func_num, Alg, runs)

% A novel cooperate coevolution framework based on difficulty and
% contribution for large-scale optimization

% parameters definition
global initial_flag
initial_flag = 0;

s1 = RandStream('mt19937ar', 'Seed', runs);
RandStream.setGlobalStream(s1);

Max_FEs = 3e6;
popsize = 100;

if strcmp(Alg, 'SHADE')
    popsize = 100;
else
    popsize = 50;
end

if benchmark == 2013
    if func_num == 13 || func_num == 14
        dim = 905;
    else
        dim = 1000;
    end
elseif benchmark == 2010
    dim = 1000;
end

[lb, ub] = getBounds(benchmark, func_num);

% load groups
[groupNum, tempFit, used_FEs] = getGroups(benchmark, func_num);

% initialization population
pop = lb + lhsdesign(popsize, dim) * (ub-lb);
fit = feval(fname, pop, func_num);
FEs = popsize + used_FEs;

% get best member and best value
[bestfit, bestid] = min(fit);
bestmem = pop(bestid, :);

fprintf('%d, %d, %d\n', func_num, size(tempFit,2));

% save the dimension index and contributions for the each groups
% OPTS is parameter of COCC
for i = 1 : groupNum
    groupSet(i).index = tempFit{i};         % dimension index
    groupSet(i).delta = 0;                  % contribution
    groupSet(i).OPTS.first = 1;
    groupSet(i).OPTS.FEs = FEs;
    groupSet(i).difficult = 1;              % 困难度
    groupSet(i).max_iter = 100;
    groupSet(i).subpop = pop(:, groupSet(i).index);     % 每个模块对应的子群
end

% parameters of optimizer
base_iter = 100;

FEs_Count = zeros(groupNum, 1);
iters = 0;

i_min = 100;
i_incre = 100;
omit = 3;

% 获得分组的初始贡献以及困难度
for i = 1 : groupNum
    [gpop, fit, used_FEs] = getSubPop(fname, func_num, groupSet(i), bestmem);
    FEs = FEs + used_FEs;
    
    [gpop, fit, bestmem, bestfit, groupSet(i), used_FEs, current_pop, current_fit, success_rate]...
      = optimizer(Alg, fname, func_num, gpop, fit, bestmem, bestfit, lb, ub, groupSet(i).max_iter, groupSet(i));
    
    FEs_Count(i) = FEs_Count(i) + groupSet(i).max_iter;
  
    groupSet(i) = evaluateDifficult(current_pop, current_fit, success_rate, groupSet(i), omit);
    FEs = FEs + used_FEs;
    iters = iters + 1;
    fprintf('%s.%2d.%2d| Gen: %d, Func Val: %f\n', Alg, func_num, runs, iters, bestfit);
end

groupIndex = groupNum;
counts = zeros(groupNum, 1);

while (FEs < Max_FEs)
    % 选择器
    [nextIndex, groupSet, order, isChange] = selector(groupIndex, groupSet, base_iter, FEs, Max_FEs, popsize, i_min, i_incre);
    
    for groupIndex = nextIndex
        dealGroup = groupSet(groupIndex);
        counts(groupIndex) = counts(groupIndex) + 1;
        % 获得子分组
        if isChange == 1
            [gpop, fit, used_FEs] = getSubPop(fname, func_num, dealGroup, bestmem);
            FEs = FEs + used_FEs;
        end

        % 计算剩余代数
        if (FEs + (dealGroup.max_iter * size(dealGroup.subpop, 1)) > Max_FEs)
            dealGroup.max_iter = ceil((Max_FEs - FEs) / popsize);
        end
        
        if FEs >= Max_FEs
            break;
        end
        
        FEs_Count(groupIndex) = FEs_Count(groupIndex) + dealGroup.max_iter;
        
        dealGroup.OPTS.FEs = FEs;
        
        % 优化器（CMA-ES，DE, SaNSDE, LSHADE）
        [gpop, fit, bestmem, bestfit, dealGroup, used_FEs, current_pop, current_fit, success_rate]...
            = optimizer(Alg, fname, func_num, gpop, fit, bestmem, bestfit, lb, ub, dealGroup.max_iter, dealGroup);
        dealGroup = evaluateDifficult(current_pop, current_fit, success_rate,dealGroup, omit);
        FEs = FEs + used_FEs;
        groupSet(groupIndex) = dealGroup;
        iters = iters + 1;
        
        fprintf('%s.%2d.%2d| Gen: %d, Func Val: %f\n', Alg, func_num, runs, iters, bestfit);

        isChange = 1;
    end
    
end

end

function [nextIdx, groupSet, max_contribution, isChange] = selector(groupIndex, groupSet, base_iter, FEs, Max_FEs, popsize, i_min, i_incre)
%% 选择所有贡献量级相同的分组
    contributions = floor(log10([groupSet.delta]));
    
    [max_contribution, ~] = max(contributions);
    nextIdx = find(contributions == max_contribution);  % 找到所有量级相同的索引
    
    min_diff = min([groupSet(nextIdx).difficult]);
    
    pre_allocate_FEs = 0;
    
    if length(nextIdx) > 1
        if min_diff == 1
            for i = nextIdx
                groupSet(i).max_iter = i_min + i_incre;
                pre_allocate_FEs = pre_allocate_FEs + groupSet(i).max_iter;
            end
        elseif min_diff == 0
            for i = nextIdx
                groupSet(i).max_iter = i_min;
                pre_allocate_FEs = pre_allocate_FEs + groupSet(i).max_iter;
            end
        else
            beta = max([1-log(min_diff/(1-min_diff)), 0]);
            for i = nextIdx
                groupSet(i).max_iter = i_min + floor(i_incre * (groupSet(i).difficult - min_diff)^beta);
                pre_allocate_FEs = pre_allocate_FEs + groupSet(i).max_iter;
            end
        end
    end
    
    % 计算剩余代数
    if (FEs + (pre_allocate_FEs * popsize) > Max_FEs)
        [~, sidx] = sort([groupSet(nextIdx).max_iter]);
        nextIdx = nextIdx(sidx);
        if nextIdx(1) == groupIndex
            isChange = 0;
        else
            isChange = 1;
        end
    else
        if ismember(groupIndex, nextIdx)
            nextIdx(nextIdx == groupIndex) = [];
            nextIdx = [groupIndex, nextIdx];
            isChange = 0;
        else
            isChange = 1;
        end
    end
    
end

function [gpop, fit, bestmem, bestfit, group, used_FEs, current_pop, current_fit, success_rate]...
      = optimizer(Alg, fname, func_num, gpop, fit, bestmem, bestfit, lb, ub, max_iter, group)
 %% CMA-ES, DE
 
dim_index  = group.index;

% best value
oldBestfit = bestfit;

[group.subpop, fit, newbestsubmem, newbestfit, group.OPTS, used_FEs, current_pop, current_fit, success_rate]...
    = feval(Alg, fname, func_num, group.subpop, gpop, fit, ub, lb, max_iter, dim_index, group.OPTS);

delta = oldBestfit - newbestfit;

if delta > 0
    bestfit = newbestfit;
    bestmem(:, dim_index) = newbestsubmem;
    group.delta = delta;
end

end

function group = evaluateDifficult(current_pop, current_fit, success_rate, group, omit)
    NP = length(current_fit);
    d = [];
    for i = 1 : NP
        fit = current_fit{i};
        subpop = current_pop{i};
        
        % 获得平均适应值
        fit_bar = mean(fit);
        [~, bestid] = min(fit);
        
        % 计算每个个体到达最优个体的距离，平均距离
        dis = pdist2(subpop, subpop(bestid, :));
        dis_bar = mean(dis);
        
        % 计算适应值-距离相关性 fitness-distance correlation (FDC)
        r = sum((fit-fit_bar).*(dis'-dis_bar)) / (sqrt(sum((fit-fit_bar).^2)) .* sqrt(sum((dis'-dis_bar).^2)));
        
        % 困难度 d = 1 - r, 第i个索引位置对应的困难度
        d = [d, abs(r)];
    end
    
    d(isnan(d) == 1) = []; %#ok<COMPNOP>
    
    if isempty(d)
        group.difficult = 1;
    else
        if omit == 1
            group.difficult = 1 - mean(d);% * success_rate;
        elseif omit == 2
            group.difficult = 1 - success_rate;
        else
            group.difficult = 1 - mean(d) * success_rate;
        end
    end
    
end

