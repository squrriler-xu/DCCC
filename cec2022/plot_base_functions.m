[X,Y] = meshgrid(-20 : 1 : 20, -20 : 1 : 20);
Z = zeros(size(X));
R = [1, 0; 0, 1];
for i = 1 : size(X, 1)
    pop = [X(i, :); Y(i, :)];
    Z(i, :) = base_functions(pop, R, 20);
end
surf(X,Y,Z)