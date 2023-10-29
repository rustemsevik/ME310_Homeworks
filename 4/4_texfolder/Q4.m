clc; clear;
% This is Q4

f = @(X, Y) (1 - X).^2 + 100*(Y - X.^2).^2;
gf = @(X, Y) [2.*X - 2 - 400.*X.*(Y - X.^2); 200.*(Y - X.^2)];

% [a, points, fvalues, iters] = SSDM(f, gf, .01, [-1.9, 2], 1000000, 1e-6, 1000);


[a, points, fvalues, iters] = SSDM(f, gf, .00076, [-1.9, 2], 1000000, 1e-6, 1000);


% [X,Y] = meshgrid(-2:0.01:2, -1:0.01:3);
% Z = (1 - X).^2 + 100*(Y - X.^2).^2;
% contour(X,Y,Z, [0 0.1 1 2 10 100 500 2500])
% hold on
% grid on
% scatter(points(:, 1), points(:, 2), "red", "filled",MarkerFaceAlpha=0.9)
% xlabel("x"); ylabel("y")
% exportgraphics(gca, "path.png", "Resolution",600)

% iterall = [,];
% pointfinall = [];
% for i = 0.00001:0.00001:0.01
%     [a, b, fvalues, iters] = SSDM(f, gf, i, [-1.9, 2], 100000, 1e-6, 1000);
%     iterall(end + 1, :) = [i, iters];
%     pointfinall(end + 1, :) = a;
% end

% plot(iterall(:,1), iterall(:,2))
% xlabel("Step size"); ylabel("Number of Iterations")
% exportgraphics(gca, "nofvsstepsize.png","Resolution",600)

% plot(fvalues, 'diamond')
% grid on
% xlim([-180 9300]); ylim([-5 280])
% ylabel("Function Value"); xlabel("Number of Iterations")
% exportgraphics(gca, "fvsiter.png","Resolution",600)


function [Pmin, points, fvals, iter] = SSDM(f, gf, h, P, maxiter, tol, divcounter)
    points = [P]; fvals = [];
    f_val = f(P(1), P(2)); 
    fvalmin = f_val; 
    Pmin = P; 
    mincounter = 0;
    
    for i = 1:maxiter
        g_vec = gf(P(1), P(2));
        g_vec_norm = g_vec / norm(g_vec);
        t_vec = g_vec_norm' * h;
        P = P - t_vec;
        points(end + 1, :) = P;
        fvals(end + 1, :) = f_val;
        f_val_old = f_val;
        f_val = f(P(1), P(2));
        
        err = abs((f_val - f_val_old) / f_val);
        if f_val < fvalmin
            fvalmin = f_val;
            Pmin = P;
            mincounter = 0;
        else
            mincounter = mincounter + 1;
        end
        if mincounter > divcounter
            break
        end

        fprintf('%d\t%.6f\t%.6f\t%.12f\t%.12f %d\n', i, P(1), P(2), f_val, err, mincounter)
        iter = i;
        if  err < tol
            break
        end
        
        Pold = P;
    end
end
