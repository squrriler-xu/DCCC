clear all; close all;

benchmark = 2013;
myfunc = 1 : 15;

%run grouping
for func=13:length(myfunc)
    global initial_flag;
    initial_flag = 0;
    %set dimension
    if benchmark == 2013
        if func == 13 || func == 14
            dim = 905;
        else
            dim = 1000;
        end
    elseif benchmark == 2010
        dim = 1000;
    end

    [lb, ub] = getBounds(benchmark, func);
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
    
    [groups, FEs] = rdg2('benchmark_func2013', func, lb, ub, dim);
    used_FEs = FEs;
    
    filename = sprintf('./result/2013/f%02d.mat',func);
    save(filename, 'groups', 'used_FEs');

end