function [groups] = getContribution(groups, cur, delta, val)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i = 1 : size(groups, 2)
    if i == cur
        continue;
    else
        groups(i).other = groups(i).other + delta;
    end
    groups(cur).val = val;
end

end

