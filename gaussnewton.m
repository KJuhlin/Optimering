function [x, funcVal, numSteps] =  gaussnewton(phi,t,y,start,tol,use_linesearch,printout,plotout)
    %setup
    
    xDim = length(start);
    rDim = length(y);
    
    phiNew = @(x) phi(x, t) - y;
    phiLinesearch = @(x) sum(phiNew(x).^2);
    
    xCurr = start;
    phiValCurr = phiNew(xCurr);
    
    crit = realmax;
    numSteps = 0;
    %Setup for the differentiation algorithm
    J = zeros(rDim, xDim);
    basis = eye(xDim);
    diffTol = 1e-6;
    
    %algorithm
    while crit > tol
        % Calculate Jacobian J and the gradient (for the stopping
        % criterion)
        for dim = 1:xDim
            J(:,dim) = (phiNew(xCurr + diffTol*basis(:,dim)) - phiNew(xCurr-diffTol*basis(:,dim)))/(2*diffTol);
        end
        
        gradient = 2*J'*phiValCurr;
        numSteps
        dir = (J'*J)\(-J'*phiValCurr);
        if (use_linesearch)
            [lambda, ~] = linesearch_armijo(phiLinesearch, xCurr, dir);
        else
            lambda = 1;
        end
        xCurr = xCurr + lambda*dir;
        
        %Check stopping criterion
        crit = norm(gradient);
        %Update step count and so on
        phiValCurr = phiNew(xCurr);
        numSteps = numSteps + 1;
        % print out
        if (printout)
            Printout()
        end
    end
    % clean up, plotout
    if (plotout)
        Plotout()
    end
    x = xCurr;
    funcVal = phiLinesearch(x);
end

