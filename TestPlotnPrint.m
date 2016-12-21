%% Tests Printout
iter = 99;
x = [1;2;3;4];
stepSize = 20;
f=10;
maxAbsRes = 3;
grad = 3;
lsIter = 4;
lambda = 1;
dGrad = 4;

Printout(iter,x,stepSize,f,maxAbsRes,grad,lsIter,lambda,dGrad)

%% Tests Plotout
[t,y] = data1;

test = @(t) 12 * exp(-3.*t);

Plotout( test,t, y)
