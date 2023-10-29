clc; clear;

f = @(x) sqrt(9.81*x/0.25) * tanh(sqrt(9.81*0.25/x)*4) - 35;
[root, iterations] = findRoot(f, [80 400], 1e-10, 1000);


function [root, iter] = findRoot(f, x0, tol, maxIter)

% Default value for maximum iteration
if nargin < 4
    maxIter = 1000;
end

% Phase 1 (searching for a bracketing interval):

s_iter = 1;
if isscalar(x0)
    fprintf("<strong> Trial       xL             xU         " + ...
        "     f(xL)           f(xU)        Method  </strong> \n")
    % initialize search width
    delta = x0/50;
    if delta == 0
        delta = 1/50;
    end
    delta_0 = delta;
    % search for bracketing interval
    sign_change = false;
    while ~sign_change % && delta < abs(x0)
        xL = x0 - delta;
        xU = x0 + delta;
        sign_change = sign(f(xL)) ~= sign(f(xU));
        fprintf('   %d\t %.6f\t %.6f\t %.6f\t %1.5e\t %s\n', ...
            s_iter, xL, xU, f(xL), f(xU), "Search")
        delta = delta_0 + 2 * delta; 
        s_iter = s_iter + 1;
        if s_iter == 1001
            break
        end
    end
    
    % check if bracketing interval was found
    if ~sign_change
        error('Unable to find bracketing interval.');
    end
else
    % use provided interval as bracketing interval
    xL = x0(1);
    xU = x0(2);
    
    % check if bracketing interval is valid
    if sign(f(xL)) == sign(f(xU))
    error('The function does not change sign in the provided interval.');
    end
end

% Phase 2 (finding the root):
fprintf("<strong> Iter\t     xL\t\t     xU\t\t    root      " + ...
    "     f(root)       Method  </strong> \n")
iter = 1;
pert = (xU - xL) / 100;
xnew = (xL + xU) / 2; % initial guess being the mid-point of the interval
while iter < maxIter
    % fist use modified secant method (MSM)
    if abs(f(xnew)) < tol
        root = xnew;
        return
    end
    x_candidate = xnew - (f(xnew) * pert) / (f(xnew + pert) - f(xnew));
    x_candidate = xnew - (f(xnew) * pert * xnew) / (f(xnew + pert * xnew) - f(xnew));
    if xL <= x_candidate && x_candidate <= xU
        xnew = x_candidate;
        fprintf('   %d\t %.6f\t %.6f\t %.6f\t %1.5e\t %s\n',...
            iter, xL, xU, xnew, f(xnew), "Modified secant")
    else % in here we have to swith to BM
        fprintf('   \t %.6f\t %.6f\t %.6f\t %1.5e\t %s\n', xL,...
            xU, x_candidate, f(x_candidate), "Modified secant (failed)")
        xnew = (xU + xL) / 2;
        fprintf('   %d\t %.6f\t %.6f\t %.6f\t %1.5e\t %s\n', iter,...
            xL, xU, xnew, f(xnew), "Bisection")
        if f(xL) * f(xnew) < 0
            xU = xnew;
            xnew = 0.5*(xL+xnew);
        elseif f(xnew) * f(xU) < 0
            xL = xnew;
            xnew = 0.5*(xnew+xU);
        end

    end
    
    % increment iterations
    iter = iter + 1;
end

% maxIter reached without finding root
error('Maximum iterations reached without finding root.');
end
