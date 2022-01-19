function [ finalgroups, used_FEs] = fast_grouping(fname, func_num, groups, LBounds, LBoundsFitVal, UBounds)
%P2P_GROUPING 此处显示有关此函数的摘要
%   此处显示详细说明
group_num = size(groups, 2);
used_FEs=0;
finalgroups={};
while group_num>0
        if group_num ==1
            finalgroups=[finalgroups,groups{1}];
            groups(1)=[];
            group_num=group_num-1;
            continue;
        end
        src_index=1;
        des_index=2:group_num;
        relation=zeros(1,group_num);
        [relation, FEs]=identify(fname, func_num, groups, src_index, des_index, relation, LBounds, LBoundsFitVal, UBounds);
        used_FEs = used_FEs + FEs;
        index=find(relation==1);
        len=length(index);
        for i=1:len
            groups{1}=[groups{1},groups{index(i)}];
        end
        groups(index)=[];
        group_num=group_num-length(index);
        if len==0 
            finalgroups=[finalgroups,groups{1}];
            groups(1)=[];
            group_num=group_num-1;
        end
end
end



