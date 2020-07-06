%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to estimate all functions of temporal Sparse Group Regression Model %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [norms,estimates] = EstimateFunctions(itrain_set,Y,S,X,options,rawEstimates)

% options
maxiter = options.maxIter;      % maximum number of iterations
tol = options.convgTol;         % stopping criterion
verbose = options.veryVerbose;  % verbose option
Y1=itrain_set.Y1; %T1
Y2=itrain_set.Y2; %T2
Y3=itrain_set.Y3; %T3
Y4=itrain_set.Y4; %T4
V=rawEstimates.V;
Yv1=Y1*V(:,1);
Yv2=Y2*V(:,2);
Yv3=Y3*V(:,3);
Yv4=Y4*V(:,4);
data=reformatY(Yv1,Yv2,Yv3,Yv4);
Y=data.Y;
% parameters
n = size(X,1);                   % number of data points
nG = size(X,2);                 % number of groups
dG = 3;                            % number of features per group
nF = dG*nG;                    % number of features

% initializations

fhat = zeros(n,nF);             % fitted values at training data points
R = Y ;                          % residuals
epsilon = 1;                       % change since last iteration
iter = 0;                             % iteration counter
obj = sum(Y.^2);                % initial objective

% fit model to training data
while (epsilon > tol && iter < maxiter)
    obj_old = obj;
    

    
    % fit additive model on partial residual
    for f = 1:nF    
        
        % group index for current function
        funInd = round(ceil(f/dG));
        
        % genotype associated with current function
        g = round(dG*(f/dG - (funInd-1)));
        
        % indicator for data points to use
        XInd = logical(X(:,funInd) == g-1);
        
        % calculate current smoother matrix
        currentS = S(XInd,XInd);
        rowSum = sum(currentS,2);
        rowSum(rowSum == 0) = 1;
        
        % update current function
        R = R + fhat(:,f).*XInd;
        partialR = R(XInd);
        fhat(XInd,f) = (currentS*partialR)./rowSum;
        R = R - fhat(:,f).*XInd;       
        
    end
    
    % prediction    
    predY =  sum(fhat,2) ;
    err = mean((predY - Y).^2);
    
    % convergence test
    obj = err;
    epsilon = abs(obj-obj_old)/obj_old;
    iter = iter + 1;

    if (verbose)
        fprintf('\tBackfitting: At iteration %d, epsilon is %f, train MSE is %f\n',iter,epsilon,err);
    end
end

% compute norms
norms = zeros(1,nG);
for j = 1:nG
    groupInds = (j-1)*dG + (1:dG);
    currFhat = fhat(:,groupInds);
    for g = 1:dG
        XInd = X(:,j) == g-1;
        if (sum(XInd) > 0)
            %norms(j) = norms(j) + sum(currFhat(XInd,g).^2)/sum(XInd); %% WEIGHTED NORM
            norms(j) = norms(j) + sum(currFhat(XInd,g).^2)/length(XInd); %% UNWEIGHTED NORM
        end
    end
    norms(j) = sqrt(norms(j));
end

% save estimates
estimates.fhat = fhat;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%