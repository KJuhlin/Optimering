function [x, funcVal, numSteps] =  gaussnewton(phi,t,y,start,tol,use_linesearch,printout,plotout)
    %setup
    
    xDim = length(start);
    rDim = length(y);
    
    phiNew = @(x) phi(x, t) - y;
    phiLinesearch = @(x) sum(phiNew(x).^2);
    
    xCurr = start;
    phiValCurr = phiNew(xCurr);
    
    crit = ones(1,5);
    numSteps = 0;
    %Setup for the differentiation algorithm
    J = zeros(rDim, xDim);
    basis = eye(xDim);
    diffTol = (1e-1)*tol;
    
    %algorithm
    while sum(abs(crit)) > tol
        % Calculate Jacobian J and the gradient (for the stopping
        % criterion)
        for dim = 1:xDim
            J(:,dim) = (phiNew(xCurr + diffTol*basis(:,dim)) - phiNew(xCurr-diffTol*basis(:,dim)))/(2*diffTol);
        end
        
        gradient = 2*J'*phiValCurr;
        
        hessianInit = 2*(J'*J);
        
        epsilon = 1;
        hessianEst = hessianInit;
        
        while cond(hessianEst) > 1000
            hessianEst = hessianInit + epsilon * eye(size(hessianInit));
            epsilon = epsilon * 4;
        end
        
        try                         % This try/catch block should never encounter errors; the check above prevents it.
            triang = chol(hessianEst);
        catch exception
            if ~strcmp(exception.identifier, 'MATLAB:posdef')
                rethrow(exception);
            end
        end
        
        rightVector = -gradient;
        d2 = triang'\rightVector;
        
        dir = triang\d2;
        if (use_linesearch)
            [lambda, lsSteps] = linesearch_armijo(phiLinesearch, xCurr, dir, tol);
        else
            lambda = 1;
            lsSteps = 0;
        end
        xCurr = xCurr + lambda*dir;
        
        %Check stopping criterion
        stepLen = norm(dir*lambda);
        crit = [stepLen crit(1:4)];
        %Update step count and so on
        phiValCurr = phiNew(xCurr);
        numSteps = numSteps + 1;
        % print out
        if (sum(isinf(phiValCurr) + isnan(phiValCurr)) > 0)
            error('Algorithm did not converge, please try different initial values or a larger tolerance.')
        end
        
        if (printout)
            %stepLen = norm(lambda*dir);
            dGrad = gradient' * dir/norm(dir);
            Printout(numSteps, xCurr, stepLen, phiLinesearch(xCurr), max(abs(phiValCurr)), norm(gradient), lsSteps, lambda, dGrad)
        end

    end
    % clean up, plotout
    if (plotout)
        plotFunc = @(t) phi(xCurr,t);
        Plotout(plotFunc, t, y);
    end
    x = xCurr;
    funcVal = phiLinesearch(x);
end

