f = @(x,y) sin(x) + 2 * cos(y);
fx = @(x,y) cos(x);
fy = @(x,y) - 2 * sin(y);
fxx = @(x,y) -sin(x);
fyy = @(x,y) - 2 * cos(y);
fxy = @(x,y) 0;

fprintf(" order          value          relative_error      true_error \n")
Taylor(f, fx, fy, fxx, fyy, fxy, 0.5, 0.51, 0.5, 0.51);

function error = Calculate_Relative_Error(xnew, xold) 
    error = abs(xnew - xold) / xnew *100;
end

function error = Calculate_Error(xnew, xold) 
    error = abs(xnew - xold);
end

function sum = Taylor(f, fx, fy, fxx, fyy, fxy, x_0, x, y_0, y)
    hx = x - x_0;
    hy = y - y_0;
    real_value = f(x, y);
    sum = 0;
    % 0th order calculation
    sum = sum + f(x_0, y_0); 
    err_rel = Calculate_Relative_Error(real_value, sum);
    err_t = Calculate_Error(real_value, sum);
    fprintf("%4d  %25.15f  %20.15e  %20.15e \n", 0, sum, err_rel, err_t)
    % 1st order calculation
    sum = sum + fx(x_0, y_0) * hx + fy(x_0, y_0) * hy;
    err_rel = Calculate_Relative_Error(real_value, sum);
    err_t = Calculate_Error(real_value, sum);
    fprintf("%4d  %25.15f  %20.15e  %20.15e \n", 1, sum, err_rel, err_t)
    % 2nd order calculation
    sum = sum + fxx(x_0, y_0) * hx^2 / ...
        2 + fxy(x_0, y_0) * hx * hy + fyy(x_0, y_0) * hy^2 / 2;
    err_rel = Calculate_Relative_Error(real_value, sum);
    err_t = Calculate_Error(real_value, sum);
    fprintf("%4d  %25.15f  %20.15e  %20.15e \n", 2, sum, err_rel, err_t)
end