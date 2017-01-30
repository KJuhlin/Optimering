function [lambda, n_iter] = linesearch_armijo(func, x, dir, tol)
    if nargin < 4
        tol = 1e-6;
    end
    epsilon = 1e-4;
    alpha = 2;
    n_iter = 1;
    lambda = 1;

    f_0 = func(x);
    df_0 = (func(x + dir*1e-8) - func(x - dir*1e-8))/(2*1e-8);
    while true
        if (func(x + dir*lambda) > f_0 + epsilon * lambda * df_0)
            lambda = lambda/2;
        elseif (func(x + dir*alpha*lambda) < f_0 + epsilon*alpha*lambda*df_0)
            lambda = 2*lambda;
        elseif (isnan(func(x + dir*lambda)) || isinf(func(x + dir*lambda)))
            lambda = tol;
            if (isnan(func(x + dir*lambda)) || isinf(func(x+dir*lambda)))
                disp('WARNING, line search: Lambda set to 0, as a small step in the chosen direction gives inf or NaN values.')
                lambda = 0;
                return;
            end
        else
            break
        end
        n_iter = n_iter + 1;
        if (n_iter > 100)
            disp('WARNING, line search: Lambda set to 0, as max number of iterations were reached.')
            lambda = 0;
            return;
        end
    end
    if isnan(func(x+lambda*dir)) || func(x+lambda*dir)>func(x)
        error('Bad job of the line search!')
    end
end

