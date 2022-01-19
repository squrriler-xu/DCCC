
% set the search range

function [XRmin, XRmax] = getBounds(benchmark, func_num)
switch benchmark
    case 2010
        if(ismember(func_num, [1, 4, 7:9, 12:14, 17:20]))
            XRmin = -100;
            XRmax = 100;
        end
        if(ismember(func_num, [2, 5, 10, 15]))
            XRmin = -5;
            XRmax = 5;
        end
        if(ismember(func_num, [3, 6, 11, 16]))
            XRmin = -32;
            XRmax = 32;
        end
    case 2013
        if(ismember(func_num, [1, 4, 7, 8, 11:15]))
            XRmin = -100;
            XRmax = 100;
        end
        if(ismember(func_num, [2, 5, 9]))
            XRmin = -5;
            XRmax = 5;
        end
        if(ismember(func_num, [3, 6, 10]))
            XRmin = -32;
            XRmax = 32;
        end
    case 2021
        if(ismember(func_num, [1, 6, 7, 8, 13, 14]))
            XRmin = -100;
            XRmax = 100;
        end
        if(ismember(func_num, [2, 3, 9, 10]))
            XRmin = -5;
            XRmax = 5;
        end
        if(ismember(func_num, [4, 5, 11, 12]))
            XRmin = -32;
            XRmax = 32;
        end
        
end

end
