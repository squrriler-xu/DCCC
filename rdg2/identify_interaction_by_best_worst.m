function [is_interactive, used_FEs, alpha, errub]=identify_interaction_by_best_worst(fname, func_num, src_group, des_group, lb, ub, lbFit)
is_interactive = 0;
used_FEs = 0;

    dim = size(lb, 2);
    ind2 = lb;
    ind2(src_group) = ub(src_group);
    val2 = feval(fname, ind2, func_num);
    ind3 = lb;
    ind3(des_group) = ub(des_group);
    val3 = feval(fname, ind3, func_num);
    ind4 = lb;
    ind4([src_group,des_group]) = ub([src_group,des_group]);
    val4 = feval(fname, ind4, func_num);
    used_FEs = used_FEs+3;
    
    muM = eps / 2;
    gamma = @(n)((n.*muM)./(1 - n.*muM));
    for i = 1 : length(lbFit)
        delta1 = lbFit(i)-val2(i);
        delta2 = val3(i)-val4(i);
        errub = gamma(dim^0.5) * max([abs(lbFit), abs(val2), abs(val3), abs(val4)]);
        alpha = abs(delta1-delta2);
        if(alpha > errub)
            is_interactive = 1;
            break;
        end
    end
end



