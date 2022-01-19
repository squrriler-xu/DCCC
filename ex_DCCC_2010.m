function ex_DCCC()
    max_run = 25;
    benchmark = 2010;
    fname = 'benchmark_func2010';
    test_func = 1 : 20;
    test = 0;
    
    mean_result = zeros(size(test_func, 2), 1);
    std_result = zeros(size(test_func, 2), 1);
    
    framework = 'DCCC';
    strategy = 'SHADE';
    
    if test == 1
        func = 1;
        result = DCCC_final(benchmark, fname, func, strategy, 3);
    else
        for func = test_func
            if ~ismember(func, test_func)
                continue;
            end
            delete(gcp('nocreate'));
            parpool('local',max_run);
            spmd(max_run)
                result = DCCC_final(benchmark, fname, func, strategy, labindex);
            end
            result = cat(1, result{1:end});
            mean_result(func) = mean(result);
            std_result(func) = std(result);
            filename=sprintf('./result/%d_%s/%s_func%02d.csv', benchmark, strategy, framework, func);
%             fil ename=sprintf('./result/2013_SaNSDE/%s_func%02d_RDG2.csv', framework, func);
            output=[mean_result(func); std_result(func); result];
            csvwrite(filename,output);
        end
        filename=sprintf('./result/%d_%s/%s.csv', benchmark, strategy, framework);
%         filename=sprintf('./result/2013_SaNSDE/%s_RDG2.csv', framework);
        output=[test_func', mean_result, std_result];
        csvwrite(filename, output);
    end
end