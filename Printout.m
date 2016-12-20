%% Loads data

f = @(x) 14*exp(- x*1.7);
[t,y] = data2;
  

gaussnewton(f,t,y,[1,1,1,1],10^-6,0,1,1)

x1=2;
x2=2;
x3=3;
x4=4;
iter=100;
step_size = 10;

%%

    fprintf('iter\tx \t\tstep \t\tsize \tf(x) \tmax(abs(r)) \tnorm(grad) \tls iters \n');
    fprintf('%2.0f \t%3.4f \t%12.4f                         \n  %12.4f \n  %12.4f \n  %12.4f\n',iter,x1,x2,x3,x4,step_size);
    
    gaussnewton(phi,t,y,start,tol,0,1,1)
    
    %%
    
%Prints the optimal parameter values (x_i i=1,2,3,4) and the total sum of
%the residuals

 f = @(x) 14*exp(-x.*1.7);

if(printout)
    % iter x step size f(x) max(abs(r)) norm(grad) ls iters lambda grad’*d/norm(d)
x1=2;
x2=2;
x3=3;
x4=4;
iter=100;
step_size = 10;


    fprintf('iter\tx \t\tstep \t\tsize \tf(x) \tmax(abs(r)) \tnorm(grad) \tls iters \n');
    fprintf('%2.0f \t%3.4f \t%12.4f                         \n  %12.4f \n  %12.4f \n  %12.4f\n',iter,x1,x2,x3,x4,step_size);
end



%Plots the function alongside the datapoints, and the residual
if(plotout)
    %Function and datapoints
    figure(1);
    plot(t,y)
    hold on
    plot(t,f(t),'r')
    hold off
    title('Function plot')
    xlabel('t')
    legend('Datapoints','F(t)')
    
    %The residual
    figure(2)
    plot(t,y-f(t))
    title('Residual plot') 
    xlabel('t')
    
end