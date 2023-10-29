val = Compute_Angle(0.25, 2, 0.02)
low = Compute_Angle(0.25 * 0.98, 2, 0.02)
up = Compute_Angle(0.25 * 1.02, 2, 0.02)

function data = Compute_Angle(a, r, percent)
    angle =  asin((1 + a) * sqrt(1 - a / (1 + a) * r^2));
    d_by_da = (-2 * a * (r^2 - 1) - r^2 + 2) / (2 * (a + 1) * ...
        sqrt((-a * r^2 + a + 1)/(a + 1))*sqrt(a*(a*(r^2 - 1) + r^2 - 2)));
    del_y = d_by_da * percent * a;
    lower_bound = angle - del_y;
    upper_bound = angle + del_y;
    data = [angle, d_by_da, del_y, lower_bound, upper_bound];
    
end