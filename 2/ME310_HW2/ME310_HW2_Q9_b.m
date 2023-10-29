clc; clear;

f1 = @(t, a) 800 - 100 * t * (cos(a) + 1);
f2 = @(t, a) 15.2 * t^2 + 80 * t * (sin(a) - 1);
f1t = @(t, a) - 100 * (cos(a) + 1);
f1a = @(t, a) 100 * t * sin(a);
f2t = @(t, a) 30.4 * t^2 + 80 * (sin(a) - 1);
f2a = @(t, a) 80 * t * cos(a);

[t, a] = NR(f1, f2, f1t, f1a, f2t, f2a, 2.5, 0.5, 1e-6, 1000)


function [x, y] = NR(f1, f2, f1x, f1y, f2x, f2y, x0, y0, eps, max_iter)

    iter = 0; error = inf;
    x = x0; y = y0;
    
    while error > eps && iter < max_iter
        J = [f1x(x, y), f1y(x, y); f2x(x, y), f2y(x, y)]; % Jacobian matrix
        F = [f1(x, y); f2(x, y)];  % function value vector
        del = J \ F;  % Newton-Raphson step
        x = x - del(1);  % update solutions
        y = y - del(2);
        error = norm(del);  % compute error
        iter = iter + 1;
    end

end
