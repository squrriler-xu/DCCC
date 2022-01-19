for func = 1 : 20
    if ismember(func, [1, 2, 3])
        groups = 1 : 1000;
        groups = num2cell(groups);
    elseif ismember(func, [4, 5, 6, 7, 8])
        if ismember(func, [4, 5, 6])
            filename=sprintf('.\DCCC\\cec2010\\datafiles\\f%02d_opm.mat', func);
        else
            filename=sprintf('.\DCCC\\cec2010\\datafiles\\f%02d_op.mat', func);
        end
        load(filename, 'p');
        groups = {p(1 : 50)}; 
        for j = 51 : 1000
            groups{j-49} = p(j);
        end
    elseif ismember(func, [9, 10, 11, 12, 13])
        if ismember(func, [9, 10, 11])
            filename=sprintf('C:\\Users\\xpl\\Documents\\MATLAB\\DifficultyAndContributionBased\\cec2010\\datafiles\\f%02d_opm.mat', func);
        else
            filename=sprintf('C:\\Users\\xpl\\Documents\\MATLAB\\DifficultyAndContributionBased\\cec2010\\datafiles\\f%02d_op.mat', func);
        end
        load(filename, 'p');
        groups = {};
        begin_idx = 1;
        end_idx = 50;
        for j = 1 : 10
            groups{j} = p(begin_idx : end_idx);
            begin_idx = end_idx + 1;
            end_idx = begin_idx + 49;
        end
        for j = 501 : 1000
            groups{j - 490} = p(j);
        end
    elseif ismember(func, [14, 15, 16, 17, 18])
        if ismember(func, [14, 15, 16])
            filename=sprintf('C:\\Users\\xpl\\Documents\\MATLAB\\DifficultyAndContributionBased\\cec2010\\datafiles\\f%02d_opm.mat', func);
        else
            filename=sprintf('C:\\Users\\xpl\\Documents\\MATLAB\\DifficultyAndContributionBased\\cec2010\\datafiles\\f%02d_op.mat', func);
        end
        load(filename, 'p');
        groups = {};
        begin_idx = 1;
        end_idx = 50;
        for j = 1 : 20
            groups{j} = p(begin_idx : end_idx);
            begin_idx = end_idx + 1;
            end_idx = begin_idx + 49;
        end
    end

end