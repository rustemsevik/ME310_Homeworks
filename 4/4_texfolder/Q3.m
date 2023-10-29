clc; clear;

% PART A

f = @(R, H) 2.*R.*H + (pi/2).*R.^2;
c = @(R, H) 2.*(H + R) + pi.*R - 300;
h = @(R) 150 - (pi / 2 + 1) .* R;

[R,H] = meshgrid(-100:1:100, -100:1:100);
Z = f(R, H);

% h = surf(R,H,Z);
% set(h,'LineStyle','none');
% print(gcf,'foo.png','-dpng','-r300');

contour(R(100:end, 100:end),H(100:end, 100:end),Z(100:end, 100:end), [6300 7000 8000])
xlabel("R"); ylabel("H");
hold on
fplot(h, [0 60], "Color","red")


% PART B

R = 0 + (60-0)*rand(1,100000);
H = h(R);
Z = f(R, H);
[M,I] = max(Z);
R(I)
H(I)

% PART c

c = @(R, H) 2.*(H + R) + pi.*R - 300;

% constraint
con = @(x)deal(c(x(1), x(2)), []);

[x, fval] = fmincon(@(x) -f(x(1), x(2)), [1, 1], [], [], [], [], [0, 0], [100, 100], con);

%results
x
