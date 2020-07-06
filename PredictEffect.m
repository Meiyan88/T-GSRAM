%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to predict snp effects on a test set using 
% the temporal Sparse Group Regression Model %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ypred,fhattest] = PredictEffect(itrain_set,Xtrain,Stest,Xtest,estimates,rawEstimates)

% parameters
ntest = size(Xtest,1);          % number of test points
nG = size(Xtest,2);             % number of groups
dG = 3;                         % number of features per group
nF = dG*nG;                     % number of features

V=rawEstimates.V;
Y1=itrain_set.Y1*V(:,1);
Y2=itrain_set.Y1*V(:,2);
Y3=itrain_set.Y1*V(:,3);
Y4=itrain_set.Y1*V(:,4);
data=reformatY(Y1,Y2,Y3,Y4);
Y=data.Y;

R = Y - sum(estimates.fhat,2);



% predict snp influence
f0=zeros(ntest,1);
fhat = zeros(ntest,nF);
for f = 1:nF    
    
    % group index for current function
    funInd = ceil(f/dG);

    % genotype associated with current function
    g = dG*(f/dG - (funInd-1));

    % indicators for data points to use
    XtrainInd = logical(Xtrain(:,funInd) == g-1);
    XtestInd = logical(Xtest(:,funInd) == g-1);
    
    % calculate current smoother matrix
    currentS = Stest(XtestInd,XtrainInd);
    rowSum = sum(currentS,2);
    rowSum(rowSum == 0) = 1;
    
    % update current function
    fhat(XtestInd,f) = currentS*(R(XtrainInd) + estimates.fhat(XtrainInd,f))./rowSum; 

end

% predict response values
if nF > 0

Ypred = sum(fhat,2);
else
    Ypred = f0 ;
end
fhattest=fhat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%