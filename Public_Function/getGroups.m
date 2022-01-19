function [group_num, grouping_result, used_FEs] = getGroups(benchmark, func_num)
% get grouping results from ideal_grouping and combine unidimensional groups

% get grouping results
groups = {};
if benchmark == 2013
%     path_grouping_result = sprintf('./ideal_grouping_2013/f%02d.mat', func_num);
    path_grouping_result = sprintf('./rdg2/result/2013/f%02d.mat', func_num);
else
    path_grouping_result = sprintf('./rdg2/result/2010/f%02d.mat', func_num);
end

load(path_grouping_result, 'groups', 'used_FEs');
% load(path_grouping_result, 'groups');
% used_FEs = 0;

% deal grouping results
single_group      = [];
grouping_result   = {};
for i = 1 : size(groups, 2)
    if(length(groups{i}) == 1)
        single_group     = [single_group, groups{i}];
    else
        grouping_result  = [grouping_result, groups(i)];
    end
end

single_num  = length(single_group);
min_size    = 50;
while(single_num > 0)
    cur_size            = min(min_size,single_num);
    grouping_result     = [grouping_result, single_group(1:cur_size)];
    single_group(1 : cur_size) = [];
    single_num          = single_num - cur_size;
end

group_num = size(grouping_result, 2);

end

