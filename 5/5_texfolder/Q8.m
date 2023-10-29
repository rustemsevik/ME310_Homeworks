clc; clear;

x = [-1, -0.96, -0.86, -0.79, 0.22, 0.5, 0.93];
y = [-1, -0.151, 0.894, 0.986, 0.895, 0.5, -0.306];

f = @(x) 0.0436499*x^6+16.0767*x^5-0.0484021*x^4-20.1005*x^3+0.0168807*x^2+5.02946*x-0.00644606;

% fplot(f, Color="blue")
% hold on
% scatter(x, y, "filled")
% xlim([-1.2 1.2]);ylim([-2, 2]); xlabel("x"); ylabel("y")
% exportgraphics(gca, "Q8.png", Resolution=600)

NDD(x, y);
Lagrng(x, y, 6);
CSI(x, y);

function table = NDD(x, y)
    n = length(x) - 1; % size of table
    
    % table
    table = zeros(n+1, n+1);
    table(:,1) = y';
    
    % Calculate the divided differences
    for j = 2:n+1
        for i = j:n+1
            table(i,j) = (table(i,j-1) - table(i-1,j-1)) / (x(i) - x(i-j+1));
        end
    end
    
end

function polynomial = Lagrng(x, y, n)
    syms X
    sum = 0;
    for i = 1:n+1
        product = y(i);
        for j = 1:n+1
            if i ~= j
                product = product * (X - x(j)) / (x(i) - x(j));
            end
        end
        sum = sum + product;
    end
    polynomial = sum;
    polynomial = vpa(simplify(polynomial));

end

function coeffs = CSI(x, y)
    n = size(x, 2) - 1;
    
    A = zeros(n);
    b = zeros(4*n, 1);
    b(1:2*n) = [y(1), repelem(y(2:end-1), 2), y(end)]';
    
    for i = 0:(n-1)
        A(2*i+1:2*i+2, 4*i+1:4*i+4) = [x(i+1)^3, x(i+1)^2, x(i+1), 1; x(i+2)^3, x(i+2)^2, x(i+2), 1];
    end
    
    for i = 0:(n-2)
        A(2*n+2*i+1:2*n+2*i+2, 4*i+1:4*i+8) = [3*x(i+2)^2, 2*x(i+2), 1, 0, -3*x(i+2)^2, -2*x(i+2), -1, 0; 6*x(i+2), 2, 0, 0, -6*x(i+2), -2, 0, 0];
    end
    
    A(end+1,1:2) = [6*x(1), 2];
    A(end+1, end-3:end) = [6*x(end), 2, 0, 0];
    coeffs = A \ b;
end