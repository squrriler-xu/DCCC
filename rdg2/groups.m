benchmark = 2013;
func_num = 8;
path_grouping_result = sprintf('./%d/f%02d.mat',benchmark, func_num);
a = unique(alpha);
b = a(~isnan(a));
j = 1;
relation = zeros(1, 1000);
tmp = [];
for i = 1 : length(b)
    if b(i) == 0
        continue;
    end
    for j = 1 : dim
        for k = i + 1 : dim
            if alpha(j, k) ~= 0
                if alpha(j, k) > errub(j, k)
                    tmp = [tmp, i, j];
                    
                end
            end
        end
        
        for tmp()
    end
end
load(path_grouping_result);