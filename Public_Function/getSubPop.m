function [gpop, fit, used_FEs] = getSubPop(fname, func_num, group, bestmem)
popsize     = size(group.subpop, 1);  % 计算子群的规模
dim_index   = group.index;

gpop        = ones(popsize, 1) * bestmem;
gpop(:, dim_index) = group.subpop;

fit = feval(fname, gpop, func_num);
used_FEs    = popsize;

end
