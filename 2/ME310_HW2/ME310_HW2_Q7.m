clc; clear;
f = @(x) log(x);
false_position(f, 0.5, 5, 0, 4)
secant(f, 0.5, 5, 0, 4)
secant(f, 5, 0.5, 0, 4)

function false_position(f, x0, x1, e, i)
    fprintf('\n*** FALSE POSITION METHOD ***\n');
    step = 1;
    condition = true;
    previous_xr = 0;
    error = 100 * e;

    while condition
        xm = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0));
        if xm ~= 0
            error = abs((previous_xr - xm) / xm * 100);
            condition = error > e;
        end
        if xm <= 0 || x0 <= 0 || x1 <= 0
            fprintf("Negative Value in ln")
            return
        end

        f_xm = f(xm);
        fprintf(['Iteration-%d, xm = %0.8f , f(xm) = %0.8f , ' ...
            'error = %0.8f\n'], step, xm, f_xm, error);

        if f(x0) * f(xm) < 0
            x1 = xm;
        elseif f(x0) * f(xm) > 0
            x0 = xm;
        else
            condition = false;
        end
        
        previous_xr = xm;
        step = step + 1;

        if step - 1 == i
            break;
        end
    end
    
    fprintf(['Last Estimate is : %0.8f, Approximate Percent Relative ' ...
        'Error : %.8f, Number of Iterations : %d , Function Value ' ...
        ': %0.8f\n'], xm, error, step-1, f_xm);
end

function secant(f, xs_old, x_s, e, i)
    fprintf('\n*** SECANT METHOD ***\n');
    step = 1;
    condition = true;
    error = 1 + e;

    while condition
        x_i = x_s - (f(x_s) * (xs_old - x_s)) / (f(xs_old) - f(x_s));
        if x_i <= 0 || x_s <= 0 || xs_old <= 0
            fprintf("Negative Value in ln")
            return
        end
        if x_i ~= 0
            error = abs((x_i - x_s) / x_i * 100);
            condition = error > e;
        end
        f_x_i = f(x_i);

        
        fprintf(['Iteration-%d, xm = %0.8f and f(xm) = %0.8f , ' ...
            'error = %0.8f\n'], step, x_i, f_x_i, error);

        xs_old = x_s;
        x_s = x_i;
        step = step + 1;
        if step -1 == i
            break;
        end
    end
    fprintf(['Last Estimate is : %0.8f, Approximate Percent Relative ' ...
        'Error : %.8f, Number of Iterations : %d , Function Value : ' ...
        '%0.8f\n'], x_s, error, step-1, f_x_i);

end
