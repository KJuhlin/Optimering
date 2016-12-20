function estimate = differentiate(phi, x, epsilon)
    if nargin < 3
        epsilon = 1e-4;
    end
    
    nDims = length(x);
    
    basis = eye(nDims);
    
    estimate = [];
    
    for dim = 1:nDims;
        estimate = [estimate (phi(x+epsilon*basis(:,dim)) - phi(x-epsilon*basis(:,dim)))/(2*epsilon)];
    end
end

