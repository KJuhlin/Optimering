function [x, funcVal, numSteps] =  gaussnewton(phi,t,y,start,tol,use_linesearch,printout,plotout,maxSteps)
    %setup
    if nargin < 9
        maxSteps = 250;
    end
    
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
        
        hessianInit = 2*(J'*J);
        
        epsilon = 0.2;
        hessianEst = hessianInit;
        
        while true
            try
                triang = chol(hessianEst);
                break;
            catch error
                if strcmp(error.identifier, 'MATLAB:posdef')
                    hessianEst = hessianInit + epsilon * eye(size(hessianInit));
                    epsilon = epsilon * 4;
                else
                    rethrow(error);
                end
            end
        end
        
        rightVector = gradient;
        d2 = triang'\rightVector;
        
        dir = triang\d2;
        if (use_linesearch)
            [lambda, lsSteps] = linesearch_armijo(phiLinesearch, xCurr, dir);
        else
            lambda = 1;
            lsSteps = 0;
        end
        xCurr = xCurr + lambda*dir;
        
        %Check stopping criterion
        stepLen = norm(dir*lambda);
        crit = stepLen/norm(xCurr);
        %Update step count and so on
        phiValCurr = phiNew(xCurr);
        numSteps = numSteps + 1;
        % print out
        if (printout)
            %stepLen = norm(lambda*dir);
            dGrad = gradient' * dir/norm(dir);
            Printout(numSteps, xCurr, stepLen, phiLinesearch(xCurr), max(abs(phiValCurr)), norm(gradient), lsSteps, lambda, dGrad)
        end
        
        %Check whether the maximum number of steps has been exceeded
        if numSteps > maxSteps
            disp(['Max number of iterations reached; Current x values are ', num2str(xCurr'), ', with functional value ', num2str(phiLinesearch(xCurr))]);
            error('Exiting function.')
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

