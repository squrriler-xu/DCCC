function fit = benchmark_func2022(x, func_num)
% The main function of the 2022 LSOPs benchmark suite
% -------------------------------- Input ----------------------------------
% x: the decision variable
% func_num: the index of the benchmark function.
% -------------------------------- Output ---------------------------------
% fit: the fitness value

global initial_flag

x = x';
if ismember(func_num, [1 : 2])
    fit = Designed_function(x, func_num);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BENCHMARK FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = Designed_function(x, func_num)
global initial_flag 
persistent xopt p s R25 R50 R100 lb ub w index_base

    [D ps] = size(x);
    if (initial_flag == 0)
        filename=sprintf('./cec2022/datafiles_2022/f%02d.mat', func_num);
        load(filename);
        initial_flag = 1;
    end

    idx = checkBounds(x, lb, ub);
    x = x-repmat(xopt,1,ps);
    fit = 0;
    ldim = 1;
    for i=1:length(s)
        if (s(i) == 25)
            f = base_functions(x(p(ldim:ldim+s(i)-1), :), R25, index_base(i));
            ldim = ldim + s(i);
        elseif (s(i) == 50)
            f = base_functions(x(p(ldim:ldim+s(i)-1), :), R50, index_base(i));
            ldim = ldim + s(i);
        elseif (s(i) == 100)
            f = base_functions(x(p(ldim:ldim+s(i)-1), :), R100, index_base(i));
            ldim = ldim + s(i);
        end
        fit = fit + w(i)*f;
%         fit = fit + f;
    end
    fit(idx) = NaN;
    if ~isempty(idx)
        warning "Some of the solutions are violating boundary constraints.";
    end
end


%------------------------------------------------------------------------------
% This function tests a given decision vector against the boundaries of a function.
%------------------------------------------------------------------------------
function indices = checkBounds(x, lb, ub)
    indices = find(sum(x > ub | x < lb) > 0);
end

