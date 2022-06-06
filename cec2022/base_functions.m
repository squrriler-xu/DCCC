function fit = base_functions(x, R, idx)

% persistent fhd

% 1. Unimodal Functions
if     (idx ==  1) fhd = str2func('f1');   % Sphere Function
elseif (idx ==  2) fhd = str2func('f2');
elseif (idx ==  3) fhd = str2func('f3');
elseif (idx ==  4) fhd = str2func('f4');
% 2. Multi-Modal Functions
%    2.1. Strong global structure:
elseif (idx ==  5) fhd = str2func('f5');
elseif (idx ==  6) fhd = str2func('f6');
elseif (idx ==  7) fhd = str2func('f7');
elseif (idx ==  8) fhd = str2func('f8');
%    2.2. Moderate or weak global structure:
elseif (idx ==  9) fhd = str2func('f9');
elseif (idx == 10) fhd = str2func('f10');
elseif (idx == 11) fhd = str2func('f11');
elseif (idx == 12) fhd = str2func('f12');
    %3. Overlapping Functions
elseif (idx == 13) fhd = str2func('f13');
elseif (idx == 14) fhd = str2func('f14');
elseif (idx == 15) fhd = str2func('f15');
    % 4. Composition Functions
elseif (idx == 16) fhd = str2func('f16');
elseif (idx == 17) fhd = str2func('f17');
elseif (idx == 18) fhd = str2func('f18');
elseif (idx == 19) fhd = str2func('f19');
elseif (idx == 20) fhd = str2func('f20');
end

fit = feval(fhd, R, x);
end

%% ------------------------- BASE FUNCTIONS ------------------------

%--------------------------------------------------------------------------
% f1 Sphere Function 
%   The simplest base function
%--------------------------------------------------------------------------
function fit = f1(R, x)
    % Done
    fit = sum(x.*x, 1);
    fit = 1e4 * fit;
end

%--------------------------------------------------------------------------
% f2 Rotated Elliptic Function
%--------------------------------------------------------------------------
function fit = f2(R, x)
    % Done
    [D ps] = size(x);
    x = R * x;
    x = T_irreg(x);
    condition = 1e+6;
    coefficients = condition .^ linspace(0, 1, D); 
    fit = 10 * coefficients * x.^2; 
end

%--------------------------------------------------------------------------
% f3 Rotated Bent Cigar Function
%--------------------------------------------------------------------------
function fit = f3(R, x)
    % Done
    [D ps] = size(x);
    x = R * x;
    x = T_asy(x, 0.5);
    fit = 1e-6 * x(1, :).^2 + sum(x(2:D, :).*x(2:D, :), 1);
    fit = 1e4 * fit;
end

%--------------------------------------------------------------------------
% f4 Rotated Discus Function
%--------------------------------------------------------------------------
function fit = f4(R, x)
    [D ps] = size(x);
    x = R * x;
    x = T_irreg(x);
    fit = 10^6 * x(1, :).^2 + sum(x(2:D, :).*x(2:D, :), 1);
    fit = 1e2 * fit;
end

%--------------------------------------------------------------------------
% f5 Different Powers Function
%--------------------------------------------------------------------------
function fit = f5(R, x)
    [D ps] = size(x);
    condition = 4;
    coefficients = 2 + condition .^ linspace(0, 1, D); 
    fit = sqrt(sum(abs(x).^repmat(coefficients', 1, ps), 1));
    fit = 1e4 * fit;
end

%------------------------------------------------------------------------------
% f6 Rotaed Rosenbrock's Function
%------------------------------------------------------------------------------
function fit = f6(R, x)
    % Done
    [D ps] = size(x);
    x = R * (2.048 * x / 100) + 1;      
    fit = sum(100.*(x(1:D-1,:).^2-x(2:D, :)).^2+(x(1:D-1, :)-1).^2, 1);
    w = log10(fit);
    w(w == -Inf) = 0;
    w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
    fit = 10.^(7 * w) .* fit;
end

%------------------------------------------------------------------------------
% f7 Rotaed Schaffers Function
%------------------------------------------------------------------------------
function fit = f7(R, x)
    [D ps] = size(x);
    y = T_diag(T_asy((R * x), 0.5), 10);
    z = sqrt(y(1:D-1, :).^2 + y(2:D, :).^2);
    fit = ((1/(D-1)) * sum(sqrt(z) + sqrt(z) .* sin(50 * z.^0.2), 1)).^2; 
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(4 * w) .* fit;
    fit = 10.^4 .* fit;
end

%------------------------------------------------------------------------------
% f8 Ackley's Function
%------------------------------------------------------------------------------
function fit = f8(R, x)
    % Done
    [D ps] = size(x);
    a = 1e10;
    b = 0.5;
    
    x = R * 0.32 * x;
    x = T_asy(x, 0.2);
%     x = T_diag(x, 2);
    fit = sum(x.^2,1);
    fit = a - a.*exp(-b.*sqrt(fit./D))-exp(sum(cos(2.*pi.*x),1)./D)+exp(1);
end

%------------------------------------------------------------------------------
% f9 Weierstrass's Function
%------------------------------------------------------------------------------
function fit = f9(R, x)
    % Done
    [D ps] = size(x);
    
    a = 0.5;
    b = 3;
    k_max = 20;
    
    x = R * (5e-3 * x);
    
    k = (0 : k_max)';
    
    temp_1 = 0;
    for i = 1 : D
        temp_1 = temp_1 + sum(a.^k .* cos(2 * pi * b.^k .* repmat((x(i, :) + 0.5), k_max+1, 1)));
    end
    temp_2 = D * sum(a.^k .* cos(pi * b.^k));
    
    fit = temp_1 - temp_2;
    fit = 10.^7 .* fit;
    
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(8 * w) .* fit;
end

%------------------------------------------------------------------------------
% f10 Rotated Griewank's Function
%------------------------------------------------------------------------------
function fit = f10(R, x)
    % Done
    [D ps] = size(x);
    
    x = R * x;
    x = T_diag(x, 100);
    i = repmat([1:D]', 1, ps);
    fit_1 = sum(x.^2/4000, 1);
    fit_2 = prod(cos(x./ sqrt(i) ), 1);
    
    fit = fit_1 - fit_2 + 1;
    
    w = log10(fit);
    w(w == -Inf) = 0;
    w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
    fit = 10.^(8 * w) .* fit;
end

%------------------------------------------------------------------------------
% f11 Rastrigin's Function, Original
%------------------------------------------------------------------------------
function fit = f11(R, x)
    % Done
    [D, ps] = size(x);
    x = 0.0512 * x;
    x = T_diag(T_asy(T_irreg(x), 0.2), 10);
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(4 * w) .* fit;
    fit = 10^6 * fit;
end

%------------------------------------------------------------------------------
% f12 Rotated Rastrigin's Function, Original
%------------------------------------------------------------------------------
function fit = f12(R, x)
    % Done
    [D, ps] = size(x);
    x = R * (0.0512 * x);
    x = T_diag(T_asy(T_irreg(x), 0.2), 10);
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(4 * w) .* fit;
    fit = 10^6 * fit;
end

%------------------------------------------------------------------------------
% f13 Non-continuous Rotated Rastrigin's Function
%------------------------------------------------------------------------------
function fit = f13(R, x)
    % Done
    [D, ps] = size(x);
    x = R * (0.0512 * x);
    x(abs(x) > 0.5) = round(2.* x(abs(x) > 0.5))/2 ;
    x = T_diag(T_asy(T_irreg(x), 0.2), 10);
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(4 * w) .* fit;
    fit = 10^6 * fit;
end

%------------------------------------------------------------------------------
% f14 Schwefel's Function
%------------------------------------------------------------------------------
function fit = f14(R, x)
    % Done
    [D, ps] = size(x);
    
    z = 10 * x;         % TODO
    
    z = z + 420.9687462275036;
    
    fit = zeros(D, ps);
    
    fit(abs(z) <= 500) = z(abs(z) <= 500).* sin(sqrt(abs(z(abs(z) <= 500))));
    fit(z > 500) = (500 - mod(z(z > 500), 500)) .* sin(sqrt(abs(500 - mod(z(z > 500), 500)))) ...
        - ((z(z > 500)-500).^2)/(10000 *D);
    fit(z < -500) = (mod(z(z < -500), 500) - 500) .* sin(sqrt(abs(mod(z(z < -500), 500) - 500))) ...
        - ((z(z < -500) + 500).^2)/(10000 *D);
    
    fit = 418.982887 * D - sum(fit, 1);
    
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(6 * w) .* fit;
    fit = 10^6 * fit;
end

%------------------------------------------------------------------------------
% f15 Rotated Schwefel's Function
%------------------------------------------------------------------------------
function fit = f15(R, x)
    % Done
    [D, ps] = size(x);
    
    z = R * 10 * x;         % TODO
    
    z = z + 420.9687462275036;
    
    fit = zeros(D, ps);
    
    fit(abs(z) <= 500) = z(abs(z) <= 500).* sin(sqrt(abs(z(abs(z) <= 500))));
    fit(z > 500) = (500 - mod(z(z > 500), 500)) .* sin(sqrt(abs(500 - mod(z(z > 500), 500)))) ...
        - ((z(z > 500)-500).^2)/(10000 *D);
    fit(z < -500) = (mod(z(z < -500), 500) - 500) .* sin(sqrt(abs(mod(z(z < -500), 500) - 500))) ...
        - ((z(z < -500) + 500).^2)/(10000 *D);
    
    fit = 418.982887 * D - sum(fit, 1);
    
%     w = log10(fit);
%     w(w == -Inf) = 0;
%     w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
%     fit = 10.^(2 * w) .* fit;
    fit = 10^5 * fit;
end

function fit = f16(R, x)
    [D, ps] = size(x);
    z = T_diag((R * 0.05 * x), 100);
    j = [1:32]';
    fit = 1;
    for i = 1 : D
        temp_z = repmat(z(i, :), 32, 1);
        a = sum(abs(2.^j .* temp_z - round(2.^j .* temp_z))./ 2.^j, 1);
        fit = fit .* ((1 + i*a) .^ (10/D^1.2));
    end
    fit = 10/D^2 * fit - 10/D^2;
    fit = 1e10 * fit;
end

%------------------------------------------------------------------------------
% f17 Lunacek Bi-Rastrigin's Function
%------------------------------------------------------------------------------
function fit = f17(R, x)
    % Done
    [D ps] = size(x);
    
    x = 0.0512 * x;

    mu0 = 2.5;
    s = 1 - (1/(2*(sqrt(D+20))-8.2));
    
    d = 1;
    mu1 = -sqrt((mu0^2-d)/s);
    
    fit = min(sum((x - mu0).^2, 1), d*D + s*sum((x - mu1).^2, 1)) + ...
      10 * (D - sum(cos(2*pi*(x - mu0)), 1)) + sum(x.^2, 1);

    fit = 10^2 * fit;
end

%------------------------------------------------------------------------------
% f18 Rotated Lunacek Bi-Rastrigin's Function
%------------------------------------------------------------------------------
function fit = f18(R, x)
    % Done
    [D ps] = size(x);
    x = T_diag(R * 0.0512 * x, 100);

    mu0 = 2.5;
    s = 1 - (1/(2*(sqrt(D+20))-8.2));
    
    d = 1;
    mu1 = -sqrt((mu0^2-d)/s);
    
    fit = min(sum((x - mu0).^2, 1), d*D + s*sum((x - mu1).^2, 1)) + ...
      10 * (D - sum(cos(2*pi*(x - mu0)), 1)) + sum(x.^2, 1);

    fit = 10^2 * fit;
end

%------------------------------------------------------------------------------
% f19 The Expanded Rosenbrock's plus Grewangk's Function
%------------------------------------------------------------------------------
function fit = f19(R, x)
    % Done
    [D ps] = size(x);
    x = 0.05 * x + 1;
    fit = 0;
    R = 1;
    for i = 1 : D-1
        temp_fit = f6(R, x(i:i+1, :));
        fit = fit + f10(R, temp_fit);
    end
    
    temp_fit = f6(R, x([D, 1], :));
    fit = fit + f10(R, temp_fit);
    
    fit = 1e8 * fit;
end

%------------------------------------------------------------------------------
% f20 The Expanded Schaffer's F6 Function
%------------------------------------------------------------------------------
function fit = f20(R, x)
    % Done
    [D ps] = size(x);
    fit = 0;
    x = x + 1;
    x = R*x;
    for i = 1 : D-1
        fit = fit + schaffer_F6(x(i, :), x(i+1, :));
    end
    fit = fit + schaffer_F6(x(1, :), x(D, :));
    
    w = log10(fit);
    w(w == -Inf) = 0;
    w = (1 ./ (1 + exp(-w)) - 0.5) * 2;
    fit = 10.^(9 * w) .* fit;
end

function g = schaffer_F6(x, y)

    g = 0.5 + ((sin(sqrt(x.^2 + y.^2))).^2 - 0.5) ./ (1+0.001*(x.^2 + y.^2)).^2;

end

%------------------------------------------------------------------------------
% This transformation function is used to break the symmetry of symmetric 
% functions.
%------------------------------------------------------------------------------
function g = T_asy(f, beta)
    [D popsize] = size(f);
    g = f;
    temp = repmat(beta * linspace(0, 1, D)', 1, popsize); 
    ind = f > 0;
    g(ind) = f(ind).^ (1 + temp(ind) .* sqrt(f(ind)));  
end


%------------------------------------------------------------------------------
% This transformation is used to create the ill-conditioning effect.
%------------------------------------------------------------------------------
function g = T_diag(f, alpha)
    [D popsize] = size(f);
    scales = repmat(sqrt(alpha) .^ linspace(0, 1, D)', 1, popsize); 
    g = scales .* f;
end


%------------------------------------------------------------------------------
% This transformation is used to create smooth local irregularities.
%------------------------------------------------------------------------------
function g = T_irreg(f)
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end

