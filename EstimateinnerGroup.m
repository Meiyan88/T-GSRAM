%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%% Function to run sub-routine of inner group temporal Sparse Group Regression Model backfitting in
% which component functions are estimated for inner group of variables %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fhat,curnorm] = EstimateinnerGroup(fhat,currGroupInds,R,S,XInd,rowSums,lambda1,lambda2,options,fhatggd)

%% parameters
% options
maxiter = options.maxIter;      % maximum number of iterations
tol = options.convgTol;         % stopping criterion

% parameters
n = size(R,1);                  % number of data points
dG = length(currGroupInds);     % size of group

%% compute curnorm
% initializations
QR = zeros(n,1);
curnorm = 0;
for g=1:dG
    rowSum = rowSums{g};
    currR=R.*XInd(:,g);
    currP = S*currR;
    P = currP(XInd(:,g))./rowSum;
    QR(XInd(:,g)) = P;
    n_g = sum(XInd(:,g));
    if (n_g > 0)
%         curnorm = curnorm + (norm(P)^2)/n_g;
curnorm = curnorm + (norm(P,'fro')^2)/n;
    end
end

 curnorm=sqrt(curnorm);
%% % estimate the component functions
if (curnorm<=lambda1) % set all functions to zero
%     fhat = zeros(n,dG);
    fhat(:,:)=eps;
else % solve equation by fixed point iteration
    iter = 0;
    epsilon = 1;
    fhat=sum(fhat,2);
%     fhat = rand(n,1);
    while (epsilon > tol && iter < maxiter)
        fhat_old = fhat;
        fhat_norm = 0;
        for g = 1:dG
            n_g = sum(XInd(:,g));
            if (n_g > 0)
%                 fhat_norm = fhat_norm + sum(fhat(XInd(:,g)).^2)/n_g; %% WEIGHTED NORM
                fhat_norm = fhat_norm + sum(fhat(XInd(:,g)).^2)/n; %% UNWEIGHTED NORM
            end
        end
        fhat_norm = sqrt(fhat_norm);
        fhat=QR./(1 + eps + lambda1/(fhat_norm+eps)+lambda2/(sqrt(fhatggd)+eps));
        epsilon = max(abs( fhat-fhat_old ));
        iter = iter + 1;
    end
    fhat=fhat./(sqrt(fhat'*fhat));
    fhatFull = zeros(n,dG);
    for g = 1:dG
        ind = XInd(:,g);
        fhatFull(ind,g) = fhat(ind);
        fhatFull(ind,g) = fhatFull(ind,g) - mean(fhat);
    end
    fhat = fhatFull; 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%