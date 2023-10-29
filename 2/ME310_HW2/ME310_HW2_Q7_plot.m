x = [1.85463498, 1.21630782, 1.05852096, 1.01616935];
y = [0.61768790, 0.19581989, 0.05687262, 0.01604002];
x = [1.85463498, 1.21630782, 0.92001335, 1.00848885];
y = [0.61768790, 0.19581989, -0.08336709, 0.00845303];

f = @(x) log(x);
hold on
grid on
fplot(f, [0.8 2])
scatter(x, y, "filled")
ylabel("f(x)"); xlabel("x")
exportgraphics(gca, "7c.png", "Resolution",300)