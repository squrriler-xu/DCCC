function [relation, used_FEs]=identify(fname, func_num, group, src_index, des_index, relation, lb, ub, lbFit)
    used_FEs=0;
    des_group=[];
    for i=1:length(des_index)
        des_group=[des_group, group{des_index(i)}];
    end
    [is_interactive, FEs]=identify_interaction_by_best_worst(fname, func_num, group{src_index}, des_group, lb, ub, lbFit);
    
    used_FEs=used_FEs+FEs;
    if is_interactive ==1
        if length(des_index)==1
               relation(des_index)=src_index;
        else
               len=length(des_index);
               mid = floor(len/2);
               [relation, FEs]=identify(fname, func_num, group, src_index, des_index(1:mid), relation, lb, ub, lbFit);
               used_FEs=used_FEs+FEs;
               [relation, FEs]=identify(fname, func_num, group, src_index, des_index(mid+1:end), relation, lb, ub, lbFit);
               used_FEs=used_FEs+FEs;
        end
    end
end
