function [groups] = getContribution(groups, cur, delta, val)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
for i = 1 : size(groups, 2)
    if i == cur
        continue;
    else
        groups(i).other = groups(i).other + delta;
    end
    groups(cur).val = val;
end

end

