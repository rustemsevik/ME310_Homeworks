clc; clear;

t = linspace(0, 10, 100);
a = linspace(-pi, pi, 100);
[t, a] = meshgrid(t, a);
func1 = 800 - 100 .* t .* (cos(a) + 1);
func2 = 15.2 .* t.^2 + 80 .* t .* (sin(a) - 1);

% h1 = surf(t, a, func2, "FaceColor",[0 0 1]);
% hold on
% h2 = surf(t, a, func1, "FaceColor",[1 0 0]);
% ylabel("\alpha"); xlabel("time"); zlabel("function values");
fdiff = func1 - func2;
C = contours(t, a, fdiff, [0 0]);
tL = C(1, 2:end);
aL = C(2, 2:end);
fL = interp2(t, a, func1, tL, aL);
% h3 = line(tL, aL, fL, 'Color', 'yellow', 'LineWidth', 3);
% legend([h1,h2,h3], {'f_{1}(t, a)', 'f_{2}(t, a)', "Line of Intersection"});
% exportgraphics(gca, "Q9a1.jpeg", "Resolution", 300)
hold off
tp = [0 10]; ap = [-pi pi]; fp = meshgrid(tp, ap);
plot3(tL, aL, fL, LineWidth=3)
hold on
surf(tp, ap, fp,'EdgeColor','none','FaceAlpha',0.5)

ylabel("\alpha"); xlabel("time"); zlabel("function values");
grid on
exportgraphics(gca, "Q9a3.jpeg", "Resolution", 300)
