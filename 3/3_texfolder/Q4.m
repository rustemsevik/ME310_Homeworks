clc; clear;
L = 0.1; Tl = 50; Tr = 100; k = 80; q = 4e+5; N = 15;

% Create tridiagonal matrix columns
f = -2 * ones(N,1);
e = ones(N,1); e(1) = 0;
g = ones(N,1); g(N) = 0;

% Create right-hand side vector
dx = 2* L / (N +1);
r = - q * dx^2 / k * ones(N,1);
r(1) = r(1) - Tl; r(N) = r(N) - Tr;

% Create a tridiagonal full matrix
A = create_array(N);

i_guess = (Tr + Tl) / 2 *ones(N, 1);

sol = linsolve(A, r)
thom = thomas_algorithm(f, e, g, r)
L = choleskyDecomposition(-A)
chol = cholesky_solution(L, -r)
seid = gauss_seidel(A, r, i_guess, 1000, 1e-7, 1)

lambda_vals = 0.1:0.01:2;

num_iters = zeros(size(lambda_vals));
for i = 1:length(lambda_vals)
    [x, iter] = gauss_seidel(A, r, i_guess, 1000, 1e-7, lambda_vals(i));
    num_iters(i) = iter;
end

plot(lambda_vals, num_iters);
xlabel('Relaxation parameter \lambda');
ylabel('Number of iterations');
exportgraphics(gca ,"lambda_vs_noi.png", 'Resolution', 300)

function x = thomas_algorithm(f, e, g, r)
    N = length(r);
    
    % Forward substitution
    for i = 2:N
        factor = e(i) / f(i-1);
        f(i) = f(i) - factor * g(i-1);
        r(i) = r(i) - factor * r(i-1);
    end
    
    % Back substitution
    x = zeros(N,1);
    x(N) = r(N) / f(N);
    for i = N-1:-1:1
        x(i) = (r(i) - g(i) * x(i+1)) / f(i);
    end

end

function A = create_array(n)
    A = -2*eye(n);
    for i = 1:n-1
        A(i,i+1) = 1;
        A(i+1,i) = 1;
    end
end

function L = choleskyDecomposition(A)
    n = size(A, 1); L = zeros(n, n);

    for k = 1:n
        L(k,k) = sqrt(A(k,k) - sum(L(k,1:k-1).^2));
        for l = k+1:n
            L(l,k) = (A(l,k) - sum(L(l,1:k-1).*L(k,1:k-1))) / L(k,k);
        end
    end
end

function x = cholesky_solution(L, b)
    n = size(L, 1); y = zeros(n, 1); x = zeros(n, 1);
    
    %Forward
    for i = 1:n
        y(i) = (b(i) - L(i,1:i-1)*y(1:i-1)) / L(i,i);
    end
    
    %Backward
    
    for i = n:-1:1
        x(i) = (y(i) - L(i+1:n,i)'*x(i+1:n)) / L(i,i);
    end

end



function [x, iter] = gauss_seidel(A, b, x ,max_iter, e, lambda)
    n = length(b);
    xi_sum = 0;
    for iter = 1:max_iter
        for i = 1:n
            for j = 1:n
                if j ~= i
                    xi_sum = xi_sum + A(i,j)*x(j);
                end
            end
            x(i) = lambda*((b(i) - xi_sum)/A(i,i)) + (1 - lambda)*x(i);
            xi_sum = 0;
        end
        residual_norm = norm(A*x - b);
        if residual_norm < e
            converged = true;
        else
            converged = false;
        end
        
        if iter ~= 1 && converged
            break;
        end
    end
    
end
