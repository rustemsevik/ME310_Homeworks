clc; clear;

g = @(c) 5 * c / (5 + 0.8 * c + c^2 + 0.2 * c^3);
gx = @(c) -(25*(2*c^3 + 5*c^2 - 25))/(c^3 + 5*c^2 + 4*c + 25)^2;
gxx = @(c) -(50*(- 3*c^5 - 15*c^4 - 21*c^3 + 150*c^2 + 375*c + 100))/(c^3 + 5*c^2 + 4*c + 25)^3;


Gold(0,10,1000, 1e-7, g)
parabolic_interpolation(g, 0, 5, 10, 1000, 1e-6, 2)
newton(g, gx, gxx, 2, 1000, 0.1)
newton_modified_derivatives(g, 0.01, 2, 1000, 0.1)


function [xopt, fx] = Gold(xl, xu, maxit, es, f)

    R = (sqrt(5) - 1) / 2; % golden ratio
    d = R * (xu - xl);
    x1 = xl + d; x2 = xu - d;
    f1 = f(x1); f2 = f(x2);
    iter = 1;
    
    if f1 > f2
        xopt = x1;
        fx = f1;
    else
        xopt = x2;
        fx = f2;
    end
    fprintf(" iter       xl   \t   xu \t\t    x1 \t\t x2 \t\t f1 \t\t f2 \t\t d\n")
    while iter < maxit
        fprintf('   %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t\n', iter, xl, xu, x1, x2, f1, f2, d)
        d = R * d;
        xint = xu - xl;
        
        if f1 > f2
            xl = x2;
            x2 = x1;
            x1 = xl + d;
            f2 = f1;
            f1 = f(x1);
        else
            xu = x1;
            x1 = x2;
            x2 = xu - d;
            f1 = f2;
            f2 = f(x2);
        end
        
        iter = iter + 1;
        
        if f1 > f2
            xopt = x1;
            fx = f1;
        else
            xopt = x2;
            fx = f2;
        end
        
        if xopt ~= 0
            ea = (1.2 - R) * abs(xint / xopt) * 100;
        end
    %     if abs(xu - xl) < 0.01
    %         break
    %     end
        if ea <= es || iter >= maxit
            break;
        end
    end

end


function [xopt, fopt] = parabolic_interpolation(f, x0, x2, x1, maxiter,tol, mode)

    f0 = f(x0); f1 = f(x1); f2 = f(x2);
    err = 100;
    x3 = 0;
    fprintf('Iter\t\tx0\t\tf0\t\t x1\t\t f1\t\t x2\t\t f2\t\t x3\t   f3\n')
    for iter = 1:maxiter
        x3old = x3;
        x3 = (f(x0)*(x1^2 - x2^2) + f(x1)*(x2^2 - x0^2) + f(x2)*(x0^2 - x1^2))/(2*f(x0)*(x1 - x2) + 2*f(x1)*(x2 - x0) + 2*f(x2)*(x0 - x1));
        f3 = f(x3);
        fprintf('%d\t %d\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', iter, x0, f0, x1, f1, x2, f2, x3, f3, err)

        if mode == 1
            x0 = x1; x1 = x2; x2 = x3;
        elseif mode == 2
            if f3 > f2 && x3 > x2 % move on with x2 x1 interval
                % Discard x0
                x0 = x2;
                x2 = x3;
                f0 = f2;
                f2 = f3;
            elseif f3 < f2 && x3 > x2 % move on with x0 x3 interval
                % Discard x1
                x1 = x3;
                f1 = f3;
            elseif f3 > f2 && x3 < x2 % move on with x0 x2 interval
                % Discard x1
                x1 = x2;
                x2 = x3;
                f1 = f2;
                f2 = f3;
            elseif f3 < f2 && x3 < x2 % move on with x0 x2 interval
                % Discard x0
                x0 = x3;
                f0 = f3;
            end
        end
        err = abs((x3 - x3old)/ x3)*100;
        
        if err < tol
            break
        end
    end
    
    % return the maximum and its location
    xopt = x3;
    fopt = f(x3);
end


function [xopt, fopt] = newton(f, fx, fxx, x, maxiter, tol)
    err = 100; 
    for i = 1:maxiter
        xnew = x - fx(x) / fxx(x);
        err = abs((x - xnew) / xnew)*100;
        fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t  %.6f\t %.6f\t\n', i, x, f(x), fx(x), fxx(x), xnew, err)
        if err < tol
            break
        end

        x = xnew;
    end
end

function [xopt, fopt] = newton_modified_derivatives(f, delta, x, maxiter, tol)
    fx = @(x) (f(x + delta*x) - f(x - delta*x)) / (2 * delta*x);
    fxx = @(x) (f(x + delta*x) - 2 * f(x) + f(x - delta*x)) / (delta*x)^2;
    err = 100;
    for i = 0:maxiter
        xnew = x - fx(x) / fxx(x);
        err = abs((x - xnew) / xnew)*100;
        fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t  %.6f\t %.6f\t\n', i, x, f(x), fx(x), fxx(x), xnew, err)
        if err < tol
            break
        end
        x = xnew;
    end
end
