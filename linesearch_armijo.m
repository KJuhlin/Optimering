function [lambda, n_iter] = linesearch_armijo(func, x, dir)
    epsilon = 1e-4;
    alpha = 2;
    n_iter = 1;
    lambda = 1;
    f_0 = func(x);

    df_0 = (func(x + dir*1e-4) - func(x - dir*1e-4))/(2*1e-4);
    while true
        if (func(x + dir*lambda) > f_0 + epsilon * lambda * df_0)
            lambda = lambda/2;
        elseif (func(x + dir*alpha*lambda) < f_0 + epsilon*alpha*lambda*df_0)
            lambda = 2*lambda;
        elseif (isnan(func(x + dir*lambda)) || isinf(func(x + dir*lambda)))
            disp('WARNING, line search: Lambda set to 0, as a unit step in the chosen direction gives inf or NaN values.')
            lambda = 0;
            return;
        else
            break
        end
        n_iter = n_iter + 1;
    end
    if isnan(func(x+lambda*dir)) || func(x+lambda*dir)>func(x)
        error('Bad job of the line search!')
    end
end

