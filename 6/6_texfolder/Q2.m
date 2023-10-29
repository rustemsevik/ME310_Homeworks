clc; clear;

f = @(x) x .* sin(x);
a = 0; b = 2*pi; 
n = 1;


Q = quad(f, a, b); % Calculate the integral using quad
result = MSTR(f, a, b, n);

while abs(result - Q) > 0.050331584227364
    n = n + 1;
    result = MSTR(f, a, b, n);
end
result
n

n = 2;
result = S13R(f, a, b, n);

while abs(result - Q) > 0.050331584227364
    n = n + 2;
    result = S13R(f, a, b, n);
end

result
n

n = 3;
result = S38R(f, a, b, n);

while abs(result - Q) > 0.050331584227364
    n = n + 3;
    result = S38R(f, a, b, n);
end

result
n

function integral = MSTR(f, a, b, n)
    % n: Number of segments
    
    h = (b - a) / n; %  width
    x = a:h:b; 
    
    integral = (h / 2) * (f(a) + f(b)); % dx area
    
    for i = 2:n
        integral = integral + h * f(x(i)); % sum dx area
    end
end

function result = S13R(f, a, b, n)
    h = (b - a) / n;
    x = a:h:b;
    y = f(x);

    result = (h / 3) * (y(1) + 4 * sum(y(2:2:n)) + 2 * sum(y(3:2:n-1)) + y(n + 1));
end

function result = S38R(f, a, b, n)
    h = (b - a) / n;
    x = a:h:b;
    y = f(x);

    result = (3 * h / 8) * (y(1) + 3 * sum(y(2:3:end-2)) + 3 * sum(y(3:3:end-1)) + 2 * sum(y(4:3:end-3)) + y(end));
end

