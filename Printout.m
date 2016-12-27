function Printout(iter,x,stepSize,f,maxAbsRes,grad,lsIter,lambda,dGrad)
%Prints some information for the optimization method

    %Local strings
    s1 = ' ';
    s2 = '  ';
    s3 = '   ';
    s4 = '    ';
    s5 = '     ';
    iterStr = [s1 'iter'];
    xStr = [s4 'x' s1];
    stepsizeStr = [s3 'step size'];
    fStr = [s2 'f(x)' s1];
    maxStr = [s2 'max(abs(r))'];
    normStr = [s2 'norm(grad)'];
    lsStr = [s2 'ls iters'];
    lambdaStr = [s2 'lambda'];
    gradStr = [s2 'grad�*d//norm(d)'];

    %Header
    fprintf('%s %s %s %s %s %s %s %s %s %s %s ', iterStr,xStr,stepsizeStr,fStr,maxStr,normStr,lsStr,lambdaStr,gradStr);

    %Vauels, only prints first x-coord
    fprintf('\n%4.0f   %4.4f    %3.4f    %4.2f     %4.4f        %4.4f     %4.0f       %4.4f      %4.4f \n', ...
            iter,x(1),stepSize,f, maxAbsRes, grad,lsIter,lambda,dGrad );
    %Prints the remaining x-coords
    for i=2:length(x);
        fprintf(' %s %3.4f \n',s5,x(i));
    end

end