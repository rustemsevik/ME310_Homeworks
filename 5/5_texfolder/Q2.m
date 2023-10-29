clc; clear

syms B

x = [1, 2, 3, 4, 5, 6, 7, 8];
y = [8.0, 12.3, 15.5, 16.8, 17.1, 15.8, 15.2, 14.0];

f = @(x) 9.77.*x.*exp(-0.215.*x);

y_fit = f(x);

St = sum((y - mean(y)).^2);
Sr = sum((y - y_fit).^2);
r2 = (St - Sr) / St

% fplot(f, [0 15], Color="blue")
% hold on
% scatter(x, y, "filled")
% ylim([0 18]); xlabel("Time[h]"); ylabel("Concentration [ng/ml]")
% exportgraphics(gca, "Q2plot.png", Resolution=600)

% Part B

% Calculate the four sums
sum1 = sum(y .* t .* exp(B * t));
sum2 = sum(t.^2 .* exp(2 * B * t));
sum3 = sum(t.^3 .* exp(2 * B * t));
sum4 = sum(y .* t.^2 .* exp(B * t));

eq = (sum1 / sum2) * sum3 - sum4;
eqa = sum1 / sum2;

B_value = vpasolve(eq == 0, B, 0)

numericA = subs(eqa, B, B_value);
A = double(numericA)

f = @(x) 9.796928.*x.*exp(-0.215087.*x);

y_fit = f(x);


St = sum((y - mean(y)).^2);
Sr = sum((y - y_fit).^2);
r2 = (St - Sr) / St

% fplot(f, [0 15], Color="blue")
% hold on
% scatter(x, y, "filled")
% ylim([0 18]); xlabel("Time[h]"); ylabel("Concentration [ng/ml]")
% exportgraphics(gca, "Q2plot2.png", Resolution=600)
