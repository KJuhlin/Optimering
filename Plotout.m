function Plotout( F,t, y)
%Plots the function alongside the datapoints, and the residual

    %Function and datapoints
    figure(1);
    plot(t,y)
    hold on
    plot(t,F(t),'r')
    hold off
    title('Function plot')
    xlabel('t')
    legend('Datapoints','F(t)')
    
    %The residual
    figure(2)
    plot(t,y-F(t))
    title('Residual plot') 
    xlabel('t')
    
end

