
% Optimize the subcomponent using SaNSDE
% The SaNSDE algorithm can be found in:
% Zhenyu Yang, Ke Tang and Xin Yao, "Self-adaptive Differential Evolution with
% Neightborhood Search", in Proceedings of the 2008 IEEE Congress on 
% Evolutionary Computation (CEC2008), Hongkong, China, 2008, pp. 1110-1116.

function [pop, fit, subbestmem, bestfit, OPTS, used_FEs, current_pop, current_fit, success_rate] = SaNSDE(fname, func_num, pop, gpop, fit, ub, lb, max_iter, dim_index, OPTS)

[popsize, dim] = size(pop);
NP = popsize;
D = dim;
% tracerst = [];

if (OPTS.first == 1)
    F = zeros(NP,1);

    linkp = 0.5;
    l1 = 1;l2 = 1;nl1 = 1;nl2 = 1;

    fp = 0.5;
    ns1 = 1; nf1 = 1; ns2 = 1; nf2 = 1;

    cc_rec = [];
    f_rec = [];

    iter = 0;

    ccm = 0.5;

    OPTS.first = 0;
else
    F      = OPTS.F;
    linkp  = OPTS.linkp;
    l1     = OPTS.l1;
    l2     = OPTS.l2;
    nl1    = OPTS.nl1;
    nl2    = OPTS.nl2;
    fp     = OPTS.fp;
    ns1    = OPTS.ns1; 
    nf1    = OPTS.nf1; 
    ns2    = OPTS.ns2; 
    nf2    = OPTS.nf2;
    cc_rec = OPTS.cc_rec;
    f_rec  = OPTS.f_rec;
    iter   = OPTS.iter;
    ccm    = OPTS.ccm;
    cc     = OPTS.cc;
end

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize DE_gbestber  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

[bestfit, ibest] = min(fit);
subbestmem = pop(ibest, :);
used_FEs = 0;

iter = 0;

current_pop = num2cell(pop, 2);
current_fit = num2cell(fit);
success_num = 0;

while iter < max_iter
    popold = pop;                   % save the old population
    
    ind = randperm(4);              % index pointer array
    
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);                
    rt = rem(rot+ind(3),NP);
    a4  = a3(rt+1);               
    rt = rem(rot+ind(4),NP);
    a5  = a4(rt+1); 
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    bm = ones(NP, 1) * subbestmem;
    
    if rem(iter,25)==0
        if (iter~=0) && (~isempty(cc_rec))
            ccm = sum(f_rec.*cc_rec)/sum(f_rec);
        end
        cc_rec = [];
        f_rec = [];
    end
    
    if rem(iter,5)==0
        cc = normrnd(ccm, 0.1, NP*3, 1);
        index = find((cc < 1) & (cc > 0));
        cc = cc(index(1:NP));
    end

    fst1 = (rand(NP,1) <= fp);
    fst2 = 1-fst1;

    fst1_index = find(fst1 ~= 0);
    fst2_index = find(fst1 == 0);

    tmp = normrnd(0.5, 0.3, NP, 1);
    F(fst1_index) = tmp(fst1_index);

    tmp = normrnd(0, 1, NP, 1) ./ normrnd(0, 1, NP, 1);
    F(fst2_index) = tmp(fst2_index);

    F = abs(F);
    
    % all random numbers < CR are 1, 0 otherwise
    aa = rand(NP,D) < repmat(cc,1,D);
    index = find(sum(aa') == 0);
    tmpsize = size(index, 2);
    for k=1:tmpsize
        bb = ceil(D*rand);
        aa(index(k), bb) = 1;
    end
        
    mui=aa;
    mpo = mui < 0.5;                % inverse mask to mui

    aaa = (rand(NP,1) <= linkp);
    aindex=find(aaa == 0);
    bindex=find(aaa ~= 0);
    
    if ~isempty(bindex)
        % mutation
        ui(bindex,:) = popold(bindex,:)+repmat(F(bindex,:),1,D).*(bm(bindex,:)-popold(bindex,:)) + repmat(F(bindex,:),1,D).*(pm1(bindex,:) - pm2(bindex,:) + pm3(bindex,:) - pm4(bindex,:));
        % crossover
        ui(bindex,:) = popold(bindex,:).*mpo(bindex,:) + ui(bindex,:).*mui(bindex,:);
    end
    if ~isempty(aindex)
        ui(aindex,:) = pm3(aindex,:) + repmat(F(aindex,:),1,D).*(pm1(aindex,:) - pm2(aindex,:));
        ui(aindex,:) = popold(aindex,:).*mpo(aindex,:) + ui(aindex,:).*mui(aindex,:);
    end
    bbb=1-aaa; 

    %-----Select which vectors are allowed to enter the new population-------
    reflect = find(ui > ub);
    ui(reflect) = ub - mod((ui(reflect) - ub), (ub - lb));
    
    reflect = find(ui < lb);
    ui(reflect) = lb + mod((lb - ui(reflect)), (ub - lb));
    
    gpop(:, dim_index) = ui;
    tempval = feval(fname, gpop, func_num);
    used_FEs = used_FEs + NP;
    a = tempval < fit;
    for i=1:NP
        if (tempval(i) <= fit(i))
            if (tempval(i) < fit(i))
                cc_rec = [cc_rec cc(i,1)];
                f_rec = [f_rec (fit(i) - tempval(i))];
            end

            pop(i,:) = ui(i,:);  
            fit(i)   = tempval(i);  
            
            l1 = l1 + aaa(i);
            l2 = l2 + bbb(i);

            ns1 = ns1 + fst1(i);
            ns2 = ns2 + fst2(i);
        else
            nl1 = nl1 + aaa(i);
            nl2 = nl2 + bbb(i); 

            nf1 = nf1 + fst1(i);
            nf2 = nf2 + fst2(i);
        end
        
            
        current_pop{i} = [current_pop{i}; ui(i, :)];
        current_fit{i} = [current_fit{i}, tempval(i)];
    end 
 
    if (rem(iter,50) == 0) && (iter~=0);
        linkp = (l1/(l1+nl1))/(l1/(l1+nl1)+l2/(l2+nl2));
        l1 = 1;l2 = 1; nl1 = 1; nl2 = 1;
        fp = (ns1 * (ns2 + nf2))/(ns2 * (ns1 + nf1) + ns1 * (ns2 + nf2));
        ns1 = 1; nf1 = 1; ns2 = 1; nf2 = 1;
    end
   
    % get best member and best value
    [temp_bestfit, ibest] = min(fit);
    temp_subbestmem = pop(ibest, :);
    
    if temp_bestfit < bestfit
        bestfit = temp_bestfit;
        subbestmem = temp_subbestmem;
        success_num = success_num + 1;
    end
    
%     tracerst = [tracerst; bestval];
    iter = iter + 1;
end

success_rate = success_num ./ max_iter;

% get best member and best value
[bestfit, bestid] = min(fit);
subbestmem = pop(bestid, :);

OPTS.F      = F;
OPTS.linkp  = linkp;
OPTS.l1     = l1;
OPTS.l2     = l2;
OPTS.nl1    = nl1;
OPTS.nl2    = nl2;
OPTS.fp     = fp;
OPTS.ns1    = ns1;
OPTS.nf1    = nf1;
OPTS.ns2    = ns2;
OPTS.nf2    = nf2;
OPTS.cc_rec = cc_rec;
OPTS.f_rec  = f_rec;
OPTS.iter   = iter;
OPTS.ccm    = ccm;
OPTS.cc     = cc;
