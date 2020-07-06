%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to select optimal groups for Ttemporal Sparse Group Regression Model %

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [groups,norms,estimates] = SelectGroups(itrain_set,Y,S,X,opts,maxgr,options,XInd_all)

%%
% clear XInd_all
% options
maxiter = options.maxIter;                  % maximum number of iterations
tol = options.convgTol;                     % stopping criterion
verbose = options.veryVerbose;              % verbose option

Y1=itrain_set.Y1; %T1
Y2=itrain_set.Y2; %T2
Y3=itrain_set.Y3; %T3
Y4=itrain_set.Y4; %T4

%% % parameters
n = size(Y,1);         % number of data points
q=size(Y,2);           %number of all times phenotype
p=size(X,2);           %number of snp
t=4;                   %time T
groupnum=opts.groupnum; %number of each group for snps
nG = length(groupnum);   % number of groups
dG = 3;           % number of features per group
nF = dG*p;        % number of features
q1 = size(Y1,2);  % number of phenotype
q2 = size(Y2,2);  % number of phenotype
q3 = size(Y3,2);  % number of phenotype
q4 = size(Y4,2);  % number of phenotype
v1 = ones(q1, 1); % initialize v1 here
v2 = ones(q2, 1); % initialize v2 here
v3 = ones(q3, 1); % initialize v3 here
v4 = ones(q4, 1); % initialize v4 here
V = [v1 v2 v3 v4]; % initialize V here
%% % initializations
f0 = zeros(n,1);                         % fitted intercept function
fhat = (zeros(n,nF));                    % fitted function values at data points
epsilon = 1;                             % change since last iteration
iter = 0;                                % iteration counter
edata=reformatY(Y1,Y2,Y3,Y4);
Y=edata.Y;
Yv1 = Y1*v1;
Yv2 = Y2*v2;
Yv3 = Y3*v3;
Yv4 = Y4*v4;
YY1 = Y1'*Y1;
YY2 = Y2'*Y2;
YY3 = Y3'*Y3;
YY4 = Y4'*Y4;
edata=reformatY(Yv1,Yv2,Yv3,Yv4);
obj = mean((edata.Y-sum(fhat,2)).^2);                        % initial objective
R=edata.Y-sum(fhat,2)-f0; % residuals
% set stopping criteria
tf = inf;
tv1 = inf;
tv2 = inf;
tv3 = inf;
tv4 = inf;

%% % generate {0,1,2} index for all SNPs
XInd_all = logical(XInd_all);
%% % precompute row sums, calculate max norm
rowSums = cell(1,nF);
maxNorm2 = 0;
for j = 1:nG
    if j>1
        groupInds = sum(groupnum(1:j-1))*dG + (1:groupnum(j)*dG); % index
    else
        groupInds =1:groupnum(j)*dG;
    end
    currNorm = 0;
    for i=1:length(groupInds)
        XInd = XInd_all(:,groupInds(i));
        if (sum(XInd) == 0)
            currS=0;
            rowSum = sum(currS,2)+eps;
        else
            currS = S(XInd,XInd);
            rowSum = sum(currS,2);
            %             rowSum(rowSum == 0) = 1;
        end
        if (sum(XInd) > 0)
            currF = currS*R(XInd)./(rowSum);
            %             currNorm = currNorm + (norm(currF)^2)/sum(XInd); %% WEIGHTED NORM
            currNorm = currNorm + (norm(currF)^2)/length(XInd); %% UNWEIGcurrNormHTED NORM
        end
        rowSums{groupInds(i)} = rowSum; 
    end
    currNorm = sqrt(currNorm);
    if (currNorm > maxNorm2)
        maxNorm2 = currNorm;
    end
end

%% % adjust lambda
lambda = opts.lambda;


%% % fit model to training data
while ((tf>tol||tv1>tol || tv2>tol || tv3>tol || tv4>tol) && iter < maxiter)

    f_old=fhat;
    
    ind = randperm(nG);
    for j = 1:nG
        if ind(j)>1
            % function indices for current group
            groupInds = sum(groupnum(1:(ind(j)-1)))*dG+(1:groupnum(ind(j))*dG);
        else
            groupInds =1:groupnum(ind(j))*dG;
        end
        % indicator for data points to use
        XInd=XInd_all(:,groupInds);
        % update functions for current group
        R = R + sum(fhat(:,groupInds).*XInd,2);
        fhat(:,groupInds) = EstimateGroup(groupInds,R,S,XInd,rowSums(groupInds),lambda,options);
        R = R - sum(fhat(:,groupInds).*XInd,2);
        
    end        
    % update v
    % -------------------------------------
    v_old1 = v1;
    v_old2 = v2;
    v_old3 = v3;
    v_old4 = v4;
    % update v
    % -------------------------------------
    di1 = updateD(v1);
    Di1 = lambda.v1*diag(di1);
    di2 = updateD(v2);
    Di2 = lambda.v1*diag(di2);
    di3 = updateD(v3);
    Di3 = lambda.v1*diag(di3);
    di4 = updateD(v4);
    Di4 = lambda.v1*diag(di4);
    dv12 = updateDV_FP(v1,v2);
    dv123 = updateDV_FP(v1,v2,v3);
    dv234 = updateDV_FP(v2,v3,v4);
    dv34 = updateDV_FP(v3,v4);
    Dv12 = lambda.v2*diag(dv12);
    Dv123 = lambda.v2*diag(dv123);
    Dv234 = lambda.v2*diag(dv234);
    Dv34 = lambda.v2*diag(dv34);
    ds = updateDs(v1,v2,v3,v4);
    DS = lambda.v3*diag(ds);
    % ------------------------
%     
    data=reformatf(fhat);
    
    f1=data.f1;
    f2=data.f2;
    f3=data.f3;
    f4=data.f4;

    % v1
    F2 = YY1+Di1+Dv12+DS;
    b2 = Y1'*sum(f1,2);
    v1 = F2\b2;
    sv1 = sqrt(v1'*YY1*v1);
    v1 = v1 ./ sv1;
    % v2
    F2 = YY2+Di2+Dv123+DS;
    b2 = Y2'*sum(f2,2);
    v2 = F2\b2;
    sv2 = sqrt(v2'*YY2*v2);
    v2 = v2 ./ sv2;
    % v3
    F2 = YY3+Di3+Dv234+DS;
    b2 = Y3'*sum(f3,2);
    v3 = F2\b2;
    sv3 = sqrt(v3'*YY3*v3);
    v3 = v3 ./ sv3;
    % v4
    F2 = YY4+Di4+Dv34+DS;
    b2 = Y4'*sum(f4,2);
    v4 = F2\b2;
    sv4 = sqrt(v4'*YY4*v4);
    v4 = v4 ./ sv4;
    V=[v1 v2 v3 v4];
    % prepare Y
    Yv1 = Y1*v1;
    Yv2 = Y2*v2;
    Yv3 = Y3*v3;
    Yv4 = Y4*v4;
    
    edata=reformatY(Yv1,Yv2,Yv3,Yv4);
    predY=sum(fhat,2);
    R=edata.Y-sum(fhat,2);
    obj_old=obj;
    norms.v2=sum(sqrt(v1.^2+v2.^2)+sqrt(v2.^2+v3.^2)+sqrt(v3.^2+v4.^2));
    norms.v3=sum(sqrt(sum(V.^2,2)));
    norms.v1=sum(sum(abs(V)));
    allNorms1 = zeros(1,p);
    for i=1:p
        group = (i-1)*dG + (1:dG);
        currFhat = fhat(:,group);
        for j=1:dG
            XInd = XInd_all(:,group(j));
            if (sum(XInd) > 0)
                allNorms1(i)=allNorms1(i)+sum(fhat(XInd,group(j)).^2)/length(XInd);
            end
        end
        allNorms1(i) = sqrt(allNorms1(i));
    end
    norms.f1=sum(allNorms1);
    allNorms1 = zeros(1,nG);
    for i=1:nG
        if j>1
            % function indices for current group
            group = sum(groupnum(1:(ind(j)-1)))*dG+(1:groupnum(ind(j))*dG);
        else
            group =1:groupnum(ind(j))*dG;
        end
        currf=fhat(:,group);
        allNorms1(i)=allNorms1(i)+sqrt(sum(sum(currf.^2))/n);
    end
    norms.f2=sum(allNorms1);
%     
    obj=mean((edata.Y-sum(fhat,2)).^2)+lambda.f1*norms.f1+lambda.f2*norms.f2+lambda.v1*norms.v1+lambda.v2*norms.v2+lambda.v3*norms.v3;
    
      
    if t > 1
        tf=abs(obj-obj_old)/obj_old;
        tv1 = max(abs(v1-v_old1));
        tv2 = max(abs(v2-v_old2));
        tv3 = max(abs(v3-v_old3));
        tv4 = max(abs(v4-v_old4));
    else
        tf = tol*10;
        tv1 = tol*10;
        tv2 = tol*10;
        tv3 = tol*10;
        tv4 = tol*10;
    end
    iter = iter+1;
end
V = [v1 v2 v3 v4];
clear norms
allNorms1 = zeros(1,p);
for i=1:p
    group = (i-1)*dG + (1:dG);
    currFhat = fhat(:,group);
    for j=1:dG
        XInd = XInd_all(:,group(j));
        if (sum(XInd) > 0)
            allNorms1(i)=allNorms1(i)+sum(fhat(XInd,group(j)).^2)/length(XInd);
        end
    end
    allNorms1(i) = sqrt(allNorms1(i));
end

% select nonzero features and groups
groups = find(allNorms1 > eps);
norms = allNorms1(groups);

% save estimates
estimates.f0 = f0;
estimates.fhat = fhat;
estimates.V=V;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%