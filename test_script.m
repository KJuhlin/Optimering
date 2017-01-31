%% Testing the linesearch

% Using the supplied testfunction
x = [0;0];
dir1 = [1;0];
dir2 = [0;1];

[lambda1, iter1] = linesearch_armijo(@test_func, x, dir1);
[lambda2, iter2] = linesearch_armijo(@test_func, x, dir2);

disp(['lambda1 = ', num2str(lambda1), ' after ', num2str(iter1), ' iterations '])
disp(['This gives the functional value ', num2str(test_func(x + lambda1*dir1))])
disp(['lambda2 = ', num2str(lambda2), ' after ', num2str(iter2), ' iterations '])
disp(['This gives the functional value ', num2str(test_func(x + lambda2*dir2))])

% With our own test function

dirPos = 1;
dirNeg = -1;
x = 0;
for a = [2 -2 5 -5 10 -10]
    func = @(x) (1-10.^a*x).^2;
    [lambdaPos, itersPos] = linesearch_armijo(func, x, dirPos);
    
    disp(['Functional value ', num2str(func(x + lambdaPos*dirPos)), ' after ', num2str(itersPos), ' iterations; pos'])
    disp(['Lambda was ', num2str(lambdaPos)])
    
    [lambdaNeg, itersNeg] = linesearch_armijo(func, x, dirNeg);
   
    disp(['Functional value ', num2str(func(x + lambdaNeg*dirNeg)), ' after ', num2str(itersNeg), ' iterations; neg'])
    disp(['Lambda was ', num2str(lambdaNeg)])
end

%% Testing the algorithm itself
% gaussnewton(function, tParams, yParams, startX, tol, linesearch, print, plot)
[t, y] = data2;
start = [5;4];
tol = 1e-4;
[x, f, steps] = gaussnewton(@phi1, t, y, start, tol, 1, 1, 1);


%% % Might want a different stopping criterion. Implementing on different branch.
[t, y] = data2;
start = [1;-2.6];
tol = 0.00001;
[x, f, steps] = gaussnewton(@phi1, t, y, start, tol, 1, 1, 1);
