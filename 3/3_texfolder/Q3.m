clc; clear;
A = [0 0.9231 0 0 0 0 0 0;
    1 -0.3846 0 0 0 0 0 0;
    0 0 0 0 -1 0 -0.8575 0;
    -1 0 -0.7809 0 0 0 0 0;
    0 0.3846 0.7809 0 1 -0.3846 0 0;
    0 -0.9231 -0.6247 0 0 0.9231 0 0;
    0 0 0.6247 1 0 0 0 0;
    0 0 0 -1 0 0 0.5145 1];
b = [1690; 3625; 0; 0; 0; 0; 0; 0];

PPGE(A, b)

function [x] = PPGE(A, b)
N = size(A,1);
x = zeros(N,1);    % Allocate memory for the solution vector

% Forward elimination
for k = 1:N-1                      % Eliminate unknown k from Eqn. k+1, k+2, … , N
    % Find the row with the largest absolute value in column k and swap with row k
    [~, maxRow] = max(abs(A(k:N, k))); % search for max in column k below row k
    maxRow = maxRow + k - 1; % adjust the row number to the original index
    if maxRow ~= k % if the row is different from the current one, swap them
        A([k, maxRow], k:N) = A([maxRow, k], k:N);
        b([k, maxRow]) = b([maxRow, k]);
    end
    for i = k+1:N                  % Modifying Eqn. i.
        factor = A(i,k) / A(k,k);
        for j = k+1:N              % Only modify {b} and the upper triangular part of [A].
                                   % Do not waste time to set the lower part of [A] to zero.
            A(i,j) = A(i,j) - A(k,j) * factor;
        end
        b(i) = b(i) - b(k) * factor;
    end
end


% Back substitution
x(N) = b(N)/A(i,i);     % Solve the last unknown using the last equation.
for i = N-1:-1:1        % This for loop runs backwards, from N-1 to 1 with increments of -1.
    sum = b(i);
    for j = i+1:N
        sum = sum - A(i,j) * x(j);
    end
    x(i) = sum / A(i,i);
end
end