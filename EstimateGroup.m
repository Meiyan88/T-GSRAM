%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to run sub-routine of temporal Sparse Group Regression Model backfitting in
% 
% which component functions are estimated for one group of variables %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fhat,curnorm] = EstimateGroup(currGroupInds,R,S,XInd,rowSums,lambda,options)

%% parameters
% options
% clear fhat
maxiter = options.maxIter;      % maximum number of iterations
tol = options.convgTol;         % stopping criterion
lambda1=lambda.f1;
lambda2=lambda.f2;
% parameters
n = size(R,1);                  % number of data points
dG = length(currGroupInds);     % size of group
sG=3;

%% compute curnorm
% initializations
QR = zeros(n,dG);
curnorm = 0;
currR=repmat(R,1,dG).*XInd;
%¼ÆËãR
for i=1:dG
    rowSum = rowSums{i};
    currP=S*currR(:,i);
    P = currP(XInd(:,i))./rowSum;
    QR(XInd(:,i),i) = P;
    n_g = sum(XInd(:,i));
    if (n_g > 0)
%                 curnorm = curnorm + (norm(P)^2)/n_g; %% WEIGHTED NORM
        curnorm = curnorm + (norm(P,'fro')^2)/n; %% UNWEIGHTED NORM
    end
end
curnorm=sqrt(curnorm);


%% % estimate the component functions
if (curnorm <= lambda2+lambda1) % set all functions to zero
    fhat = zeros(n,dG);
else % solve equation by fixed point iteration
    iter = 0;
    epsilon = 1;
    fhat=rand(n,dG).*XInd;
    while (epsilon > tol && iter < maxiter)
        fhat_old = fhat;
        fhatggd=0;
        for g = 1:dG
            n_g = sum(XInd(:,g));
            if (n_g > 0)
%                 fhatggd = fhatggd + sum(fhat(XInd(:,g)).^2)/n_g; %% WEIGHTED NORM
                fhatggd = fhatggd + sum(fhat(XInd(:,g)).^2)/n; %% UNWEIGHTED NORM
            end
        end

        for i=1:dG/sG
            groupInds=(i-1)*sG+(1:sG);
            Xind=XInd(:,groupInds);
            R = R + sum(fhat(:,groupInds).*Xind,2);
            fhat(:,groupInds)= EstimateinnerGroup(fhat(:,groupInds),groupInds,R,S,Xind,rowSums(groupInds),lambda1,lambda2,options,fhatggd);
            R = R - sum(fhat(:,groupInds).*Xind,2);
        end
        epsilon = max(max(abs(fhat_old - fhat)));
        iter = iter + 1;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%