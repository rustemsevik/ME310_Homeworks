% METU ME310 Spring 2023
% Author: Muhammed Rustem Sevik 2446888

% This code finds the first n postive roots of a given function using
% nested incremental search

% Note that this code contains two functions one of them is modified 
% version of the code on OdtuClass and other one is to make the original 
% code perecise

f = @(x) sin(x.^2) + exp(-0.5*x);

IncrementalSearch(f, 0, 100, 0.01, 5, 3, 0.1)

function [roots] = Incremental(f, xMin, xMax, dx, nRoot)

foundRoots = 0;  % Counter for roots that are found
x = xMin;

while foundRoots < nRoot && x <= xMax
    f1 = f(x);
    f2 = f(x+dx);
    
    if f1*f2 < 0  % Sign change indicates that there is a root in [x,x+dx]
        foundRoots = foundRoots + 1;

        % To locate the root in [x,x+dx] we can perform linear
        % interpolation between points (x,f1) and (x+dx,f2).
        roots(foundRoots,1:2) = [x x+dx];
    end
    x = x + dx;
end

end  % of the function


function [final_roots] = IncrementalSearch(f, xMin, xMax, dx, nRoot, ...
    nSearch, fraction)
    roots = Incremental(f, xMin, xMax, dx, nRoot);
    a = size(roots);
    final_roots = zeros(a(1,1),1);

    for i = 1 : a(1,1)
        dx_local = dx;
        
        for j = 1 : (nSearch - 1)

            dx_local = dx_local * fraction;
            roots(i,1:2) = Incremental(f,roots(i,1),roots(i,2),dx_local,1);
        
        end % of j for loop

        final_roots(i,1) = (abs(f(roots(i,1)))*(roots(i,2))+ ...
            abs(f(roots(i,2)))*roots(i,1)) / ...
            (abs(f(roots(i,1)))+abs(f(roots(i,2))));
    
    end % of i for loop

end % of the function














