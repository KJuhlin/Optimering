% Setup

clear all
close all
clc

ofile = fopen('testresults/res_large.txt', 'w');

[t1, y1] = data1;
[t2, y2] = data2;

tol = 1e-4;

% phi1: x1 * exp(-x2 * t)

x1_inits = 1:10;
x2_inits = 1:5;
x3_inits = 1:10;
x4_inits = 1:5;

use_linesearch = 0;

data1_phi1_limit = 10;
data2_phi1_limit = 250;
data1_phi2_limit = 1;
data2_phi2_limit = 10;


% phi2: x1*exp(-x2*t) + x3*exp(-x4*t)

% Testing data set 1 without line search

% First function

convergent = 0;
total = 0;
linesearch_fails = 0;

fprintf(ofile, '%s\n', 'Data set 1, no line search, phi1');
fprintf(ofile, '%10s %10s %12s %8s %8s %10s\n', 'x1_init', 'x2_init', 'Iterations', 'x1:', 'x2:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        x_init = [x1; x2];
        try
            [x_res, f_val, n_steps] = gaussnewton(@phi1, t1, y1, x_init, tol, use_linesearch, 0, 0);
            fprintf(ofile, '%10d %10d %12d %8.5f %8.5f %10.5f\n', x1, x2, n_steps, x_res(1), x_res(2), f_val);
        catch exception
            fprintf(ofile, 'No convergence\n');
        end
        if (f_val < data1_phi1_limit)
            convergent = convergent + 1;
        end
        total = total + 1;
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

% Second function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 1, no line search, phi2');
fprintf(ofile, '%10s %10s %10s %10s %12s %8s %8s %8s %8s %10s\n', 'x1_init', 'x2_init', 'x3_init', 'x4_init', ...
    'Iterations', 'x1:', 'x2:', 'x3:', 'x4:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        for x3 = x3_inits
            for x4 = x4_inits
                x_init = [x1; x2; x3; x4];
                try
                    [x_res, f_val, n_steps] = gaussnewton(@phi2, t1, y1, x_init, tol, use_linesearch, 0, 0);
                    fprintf(ofile, '%10d %10d %10d %10d %12d %8.5f %8.5f %8.5f %8.5f %10.5f\n', ...
                        x1, x2, x3, x4, n_steps, x_res(1), x_res(2), x_res(3), x_res(4), f_val);
                catch exception
                    fprintf(ofile, 'No convergence\n');
                end
                if (f_val < data1_phi2_limit)
                    convergent = convergent + 1;
                end
                total = total + 1;
            end
        end
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

% Testing data set 1 with line search

use_linesearch = 1;

% First function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 1, with line search, phi1');
fprintf(ofile, '%10s %10s %12s %8s %8s %10s\n', 'x1_init', 'x2_init', 'Iterations', 'x1:', 'x2:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        x_init = [x1; x2];
        try
            [x_res, f_val, n_steps] = gaussnewton(@phi1, t1, y1, x_init, tol, use_linesearch, 0, 0);
            fprintf(ofile, '%10d %10d %12d %8.5f %8.5f %10.5f\n', x1, x2, n_steps, x_res(1), x_res(2), f_val);
        catch exception
            fprintf(ofile, 'No convergence\n');
        end
        if (f_val < data1_phi1_limit)
            convergent = convergent + 1;
        end
        total = total + 1;
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

% Second function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 1, with line search, phi2');
fprintf(ofile, '%10s %10s %10s %10s %12s %8s %8s %8s %8s %10s\n', ...
    'x1_init', 'x2_init', 'x3_init', 'x4_init', 'Iterations', 'x1:', 'x2:', 'x3:', 'x4:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        for x3 = x3_inits
            for x4 = x4_inits
                x_init = [x1; x2; x3; x4];
                try
                    [x_res, f_val, n_steps] = gaussnewton(@phi2, t1, y1, x_init, tol, use_linesearch, 0, 0);
                    fprintf(ofile, '%10d %10d %10d %10d %12d %8.5f %8.5f %8.5f %8.5f %10.5f\n', ...
                            x1, x2, x3, x4, n_steps, x_res(1), x_res(2), x_res(3), x_res(4), f_val);
    
                catch err
                    fprintf(ofile, 'No convergence\n');
                end
                 if (f_val < data1_phi2_limit)
                    convergent = convergent + 1;
                end
                total = total + 1;
            end
        end
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);


% Testing data set 2 without line search

use_linesearch = 0;

% First function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 2, no line search, phi1');
fprintf(ofile, '%10s %10s %12s %8s %8s %10s\n', ...
    'x1_init', 'x2_init', 'Iterations', 'x1:', 'x2:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        x_init = [x1; x2];
        try
            [x_res, f_val, n_steps] = gaussnewton(@phi1, t2, y2, x_init, tol, use_linesearch, 0, 0);
            fprintf(ofile, '%10d %10d %12d %8.5f %8.5f %10.5f\n', x1, x2, n_steps, x_res(1), x_res(2), f_val);
        catch exception
            fprintf(ofile, 'No convergence\n');
        end
        if (f_val < data2_phi1_limit)
            convergent = convergent + 1;
        end
        total = total + 1;
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

% Second function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 2, no line search, phi2');
fprintf(ofile, '%10s %10s %10s %10s %12s %8s %8s %8s %8s %10s\n', ...
    'x1_init', 'x2_init', 'x3_init', 'x4_init', 'Iterations', 'x1:', 'x2:', 'x3:', 'x4:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        for x3 = x3_inits
            for x4 = x4_inits
                x_init = [x1; x2; x3; x4];
                try
                    [x_res, f_val, n_steps] = gaussnewton(@phi2, t2, y2, x_init, tol, use_linesearch, 0, 0);
                    fprintf(ofile, '%10d %10d %10d %10d %12d %8.5f %8.5f %8.5f %8.5f %10.5f\n', ...
                        x1, x2, x3, x4, n_steps, x_res(1), x_res(2), x_res(3), x_res(4), f_val);
                catch exception
                    fprintf(ofile, 'No convergence\n');
                end
                if (f_val < data2_phi2_limit)
                    convergent = convergent + 1;
                end
                total = total + 1;
            end
        end
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);


% Testing data set 2 with line search
use_linesearch = 1;

% First function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 2, with line search, phi1');
fprintf(ofile, '%10s %10s %12s %8s %8s %10s\n', 'x1_init', 'x2_init', 'Iterations', 'x1:', 'x2:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        x_init = [x1; x2];
        try
            [x_res, f_val, n_steps] = gaussnewton(@phi1, t2, y2, x_init, tol, use_linesearch, 0, 0);
            fprintf(ofile, '%10d %10d %12d %8.5f %8.5f %10.5f\n', x1, x2, n_steps, x_res(1), x_res(2), f_val);
        catch exception
            fprintf(ofile, 'No convergence\n');
        end
        if (f_val < data2_phi1_limit)
            convergent = convergent + 1;
        end
        total = total + 1;
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

% Second function

convergent = 0;
total = 0;

fprintf(ofile, '%s\n', 'Data set 2, with line search, phi2');
fprintf(ofile, '%10s %10s %10s %10s %12s %8s %8s %8s %8s %10s\n', ...
    'x1_init', 'x2_init', 'x3_init', 'x4_init', 'Iterations', 'x1:', 'x2:', 'x3:', 'x4:', 'Func val');

for x1 = x1_inits
    for x2 = x2_inits
        for x3 = x3_inits
            for x4 = x4_inits
                x_init = [x1; x2; x3; x4];
                try
                    [x_res, f_val, n_steps] = gaussnewton(@phi2, t2, y2, x_init, tol, use_linesearch, 0, 0);
                    fprintf(ofile, '%10d %10d %10d %10d %12d %8.5f %8.5f %8.5f %8.5f %10.5f\n', ...
                        x1, x2, x3, x4, n_steps, x_res(1), x_res(2), x_res(3), x_res(4), f_val);
                catch exception
                    fprintf(ofile, 'No convergence\n');
                end

                if (f_val < data2_phi2_limit)
                    convergent = convergent + 1;
                end
                total = total + 1;
            end
        end
    end
end

fprintf(ofile, '%d out of %d runs converged\n\n', convergent, total);

fprintf(ofile, '%d linesearch failures encountered in 10200 runs\n', linesearch_fails);

fclose(ofile);

% To print:
% What run
% Number of iterations
% x-values
% Functional value