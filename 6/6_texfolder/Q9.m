clc; clear;

t = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48];
r = 1000*[18.80, 18.86, 18.95, 19.04, 19.15, 19.26, 19.38, 19.50, 19.62, 19.74, 19.87, 19.99, 20.12];
theta = [0.785, 0.779, 0.770, 0.759, 0.748, 0.735, 0.722, 0.707, 0.693, 0.677, 0.661, 0.645, 0.628];
h = 4; l = length(t);

v = @(dr, dt, r, t) sqrt(dr.^2 + (r.*dt).^2);
a = @(dr, ddr, dt, ddt, r, t) sqrt((ddr - r.*dt.^2).^2 + (r.*ddt + 2.*dr.*dt).^2);

dr = 0; dt = 0; ddr = 0; ddt = 0;
vs = zeros(l, 1); 
as = zeros(l, 1);

for i = 1:l
    if i == 1 
        dr = (-r(i+2) + 4*r(i+1) - 3*r(i)) / (2*h);
        dt = (-theta(i+2) + 4*theta(i+1) - 3*theta(i)) / (2*h);
        ddr = (-r(i+3) + 4*r(i+2) - 5*r(i+1) + 2*r(i)) / h^2;
        ddt = (-theta(i+3) + 4*theta(i+2) - 5*theta(i+1)+2*theta(i)) / h^2;
    elseif i == l
        dr = -(-r(i-2) + 4*r(i-1) - 3*r(i)) / (2*h);
        dt = -(-theta(i-2) + 4*theta(i-1) - 3*theta(i)) / (2*h);
        ddr = -(-r(i-3) + 4*r(i-2) - 5*r(i-1) + 2*r(i)) / h^2;
        ddt = -(-theta(i-3) + 4*theta(i-2) - 5*theta(i-1)+2*theta(i)) / h^2;
    else
        dr = (r(i+1) - r(i-1)) / (2*h);
        dt = (theta(i+1) - theta(i-1)) / (2*h);
        ddr = (r(i+1) -2*r(i) + r(i-1))/h^2;
        ddt = (theta(i+1) -2*theta(i) + theta(i-1))/h^2;
    end
    vs(i) = v(dr, dt, r(i), t(i));
    as(i) = a(dr, ddr, dt, ddt, r(i), t(i));
end

figure;
plot(t, vs); xlabel("t [s]"); ylabel("velocity [m/s]");
axis('padded')
% exportgraphics(gcf, "9_1.eps")

figure;
plot(t, as); xlabel("t [s]"); ylabel("acceleration [m/s^{2}]");
axis('padded')
% exportgraphics(gcf, "9_2.eps")

