function generateDatafiles
% Generate configuration data for each function

D = 1000;               % the dimension of the benchmark functions
r = [0, 0.05, 0.05, 0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0];
                        % the coupling degree of the benchmark functions
type = [1, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 2, 2, 2, 3];
                        % the type of coupling topology of the benchmark functions
k = 1; mu = [0.1, 0.2, 0.3, 0.4, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0, 0.4, 0.5, 1.0, 1.2];
                        % the parameters of the weight

num_benchmark = 2;     % the number of benchmark functions
num_base = 20;          % the number of base functions
num_components = 20;    % the number of the components of the benchmark functions

for func = 3
    % Fixed the random number seed
    s1 = RandStream('mt19937ar', 'Seed', func);
    RandStream.setGlobalStream(s1);
    
    % Generate the index of the base functions used by the current benchmark function
    index_base = selectBaseFunction(num_components, num_base);
    
    % Generate the sets of decision variables for each component, a set of private variables and shared variables.
    [I_all, S] = allocateVariable(D, r(func), num_components, type(func));
    
    % Generate random weights
    w = zeros(1, num_components);
    for i = 1 : num_components
        w(i) = 10 ^ (2 * normrnd(0, 1));
    end
    
    % Generate random the type of coupling: conforting and conflicting
    
    if func == 1
        b = ones(1, num_components);
    else
        b = rand(1, num_components);
        b(b < 0.5) = 0;
        b(b >= 0.5) = 1;
    end
    
    % The upper and lower boundary
    lb = -100;
    ub = 100;
    
    % Generate random shift vector
    o = cell(num_components - sum(b)+1, 1);
    xopt = cell(num_components, 1);
    
    comf_set = [];
    bidx = 1;
    for i = 1 : num_components
        if b(i) == 1
            comf_set = [comf_set, I_all{i}];
        else
            o{bidx} = computeShiftVector(S(i), lb, ub);
            xopt{i} = o{bidx};
            bidx = bidx + 1;
        end
    end
    comf_set = unique(comf_set, 'stable');
    o{bidx} = computeShiftVector(length(comf_set), lb, ub);
    
    for i = 1 : num_components
        if b(i) == 1
            for j = I_all{i}
                xopt{i} = [xopt{i}; o{bidx}(comf_set == j)];
            end
        end
    end
    xopt_local = zeros(D, 2);
    for i = 1: 2
        xopt_local(:, i) = computeShiftVector(D, lb, ub);
    end    
    % Generate random rotation matrix
    R = cell(num_components, 1);
    for i = 1 : num_components
        R{i} = computeRotation(S(i));
    end
    
%     groups = I_all';
%     filename=sprintf('./ideal_grouping/f%02d.mat', func);
%     save(filename, 'groups');
    
    filename=sprintf('./datafiles_2022/f%02d.mat', func);
    save(filename, 'index_base', 'I_all', 'S', 'w', 'b', 'lb', 'ub', 'o', 'R', 'xopt', 'xopt_local');
end

for func = [12, 13]
    s1 = RandStream('mt19937ar', 'Seed', func);
    RandStream.setGlobalStream(s1);
    
    % The upper and lower boundary
    lb = -100;
    ub = 100;
    o = computeShiftVector(D, lb, ub);
    
    groups = cell(1, 20);
    ldim = 1;
    for i = 1 : 20
        groups{i} = ldim : ldim+49;
        ldim = ldim + 50;
    end
%     filename=sprintf('./ideal_grouping/f%02d.mat', func);
%     save(filename, 'groups');
    
    filename=sprintf('./datafiles_2022/f%02d.mat', func);
    save(filename, 'lb', 'ub', 'o');
end

end

function [index_base] = selectBaseFunction(num_components, num_base)
% Generate a random index of the base function corresponding to each
% component in the benchmark function.
% -------------------------------- Input ----------------------------------
% num_components: the number of components contained in the benchmark 
% function
% num_base: the number of base functions
% -------------------------------- Output ---------------------------------
% index_base: the index of the base functions used by the benchmark
% funcions

max_num = 2;        % Upper limit of reuse of a base function
index_base = randi(num_base, 1, num_components);

for i = 1 : num_components
    while sum(index_base == index_base(i)) > max_num
        index_base(i) = randi(num_base);
    end
end

end

function [I_all, S] = allocateVariable(D, r, num_components, type)
% Randomly allocate decision variables to each component.
% -------------------------------- Input ----------------------------------
% D: the dimension of the benchmark function
% r: the coupling degree of the benchmark function
% num_components: the number of components contained in the benchmark 
% function
% type: the type of coupling topology (1: sequential, 2: ring, 3: random)
% -------------------------------- Output ---------------------------------
% I_all: record the decision variables owned by each component
% I_priv: record the prviate decision variables owned by each component
% I_shar: record the shared decision variables owned by each component
% S: record the size of each component

num_priv = floor((1 - r) * D);     % the number of private variables
num_shar = 1000 - num_priv;           % the number of shared variables

I_all_set  = randperm(D);   % shuffle the index of decision variable
I_priv_set = I_all_set(1 : num_priv);   % the set of private variables
I_shar_set = I_all_set(num_priv+1 : end);% the set of shared variables

I_all  = cell(num_components, 1);
I_priv = cell(num_components, 1);
I_shar = cell(num_components, 1);

S_priv = zeros(num_components, 1);
S_shar = zeros(num_components, 1);

if num_priv ~= 0
    index_divide = randperm(num_priv, num_components-1);
    index_divide = sort(index_divide);
    
    ldim = 1;
    for i = 1 : num_components-1
        I_priv{i} = I_priv_set(ldim : index_divide(i));
        ldim = index_divide(i)+1;
        S_priv(i) = length(I_priv{i});
    end
    I_priv{num_components} = I_priv_set(ldim : end);
    S_priv(num_components) = length(I_priv{num_components});
    
    % lightly balance the size of different components. If the size of the component is less than 25, 
    % then the component is balanced by taking some variables from the largest component to 25
    if r < 0.2        
        for i = 1 : num_components
            min_size = 25 + randi(25);
            if S_priv(i) < min_size
                num_lack = min_size - S_priv(i);
                [~, j] = max(S_priv);
                I_priv{i} = [I_priv{i}, I_priv{j}(1 : num_lack)];
                S_priv(i) = min_size;
                I_priv{j} = I_priv{j}(num_lack+1 : end);
                S_priv(j) = S_priv(j) - num_lack;
            end
        end
    end
end

if num_shar ~= 0
    B = zeros(num_shar, num_components);
    if type == 1        % sequential coupling topology
        j = 1;
        for i = 1 : num_shar
            B(i, j) = 1;
            if j == 1
                B(i, j+1) = 1;
            elseif j == num_components
                B(i, j-1) = 1;
            else
                if rand < 0.5
                    B(i, j+1) = 1;
                else
                    B(i, j-1) = 1;
                end
            end
            
            if j == num_components
                j = 1;
            else
                j = j + 1;
            end
        end
    elseif type == 2    % ring coupling topology
        j = 1;
        for i = 1 : num_shar
            B(i, j) = 1;
            if j == 1
                if rand < 0.5
                    B(i, j+1) = 1;
                else
                    B(i, num_components) = 1;
                end
            elseif j == num_components
                if rand < 0.5
                    B(i, 1) = 1;
                else
                    B(i, j-1) = 1;
                end
            else
                if rand < 0.5
                    B(i, j+1) = 1;
                else
                    B(i, j-1) = 1;
                end
            end
            
            if j == num_components
                j = 1;
            else
                j = j + 1;
            end
        end
    elseif type == 3    % random coupling topology
        rand_matrix = rand(num_shar, num_components);
        B(rand_matrix < 0.1) = 1;
        
        for i = 1 : num_shar    % ensure that every variable is a shared variable
            if sum(B(i, :)) < 2
                B(i, :) = 0;
                B(i, randperm(num_components, 2)) = 1;
            end
        end
        
        % each component of the fully coupled problem is ensured to contain at least 25 variables
        index_lack = find(sum(B, 1) < 25);
        for i = 1 : length(index_lack)
            num_supply = 25 - (sum(B(:, index_lack)) + S_priv(index_lack));
            if num_supply > 0
                B(randperm(num_shar, num_supply), index_empty(i)) = 1;
            end
        end
        
    end
    
    for i = 1 : num_components
        I_shar{i} = I_shar_set(B(:, i) == 1);
        S_shar(i) = length(I_shar{i});
    end
end

S = S_priv + S_shar;

for i = 1 : num_components
    I_all{i} = [I_priv{i}, I_shar{i}];
end

end

function s = computeShiftVector(D, lb, ub)
% Randomly generate a shift vector to the component.
    s = zeros(D, 1);
    hw = (ub - lb) / 2.0;
    middle = lb + hw;
    for i=1:D
        s(i) = middle + randn * hw;
        while((s(i) < 0.8 * lb) || (s(i) > 0.8 * ub))
            s(i) = middle + randn * hw;
        end
    end
end


function [Q R] = computeRotation(D)
% Randomly generate a rotation matrix to the component.
    A = randn(D);
    Q = zeros(D);
    R = zeros(D);
    for j = 1:D
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end


