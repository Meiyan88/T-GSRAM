%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performs temporal Sparse Group Regression Model  %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data %
%%
clc
clear
addpath('./TSGRAM/');
dataDir='D:\cxm\tulip\TY\dataDir';
outputDir='D:\cxm\tulip\TY\outputDir\';%每次循环都要修改
load options.mat
load data.mat
% clear snps
opts.groupidx = group_idx;
opts.groupnum = group_num;
%get normalization
Y1=getNormalization(Y1);
Y2=getNormalization(Y2);
Y3=getNormalization(Y3);
Y4=getNormalization(Y4);
n=size(Y1,1);
numcv=5;%10
ind=crossvalind('Kfold',n,numcv);
lambdaVals1 = sort(options.lambdaVals.v,'descend');

disp('Begin cross validition...');
disp('===========================');
%
%%  run algothm
for k=1:5
    testind=ind==k;
    trainind=~testind;
    train.Y1=Y1(trainind,:);train.Y2=Y2(trainind,:);train.Y3=Y3(trainind,:);train.Y4=Y4(trainind,:);
    test.Y1=Y1(testind,:);test.Y2=Y2(testind,:);test.Y3=Y3(testind,:);test.Y4=Y4(testind,:);
    agesTrain=ages(trainind,:);snpsTrain=SNPs(trainind,:);
    agesTest=ages(testind,:);snpsTest=SNPs(testind,:);
    for j1=1:length(lambdaVals1)
        opts.lambda.f1=lambdaVals1(j1); % G21 norm
        for j2=1:length(lambdaVals1)
            opts.lambda.f2=lambdaVals1(j2);%L21 norm, individual sparsity
            for j3=1:length(lambdaVals1)
                opts.lambda.v1=lambdaVals1(j3);% L1-norm, individual sparsity
                for j4=1:length(lambdaVals1)
                    opts.lambda.v2=lambdaVals1(j4);%F21-norm,time-consistent norm
                    for j5=1:length(lambdaVals1)
                        opts.lambda.v3=lambdaVals1(j5);%  L21-norm, individual across tasks
                        ind1=crossvalind('Kfold',n,numcv);
                        for i=1:numcv
                            itestind=ind1==i;
                            itrainind=~testind;
                            itrain.Y1=Y1(trainind,:);itrain.Y2=Y2(trainind,:);itrain.Y3=Y3(trainind,:);itrain.Y4=Y4(trainind,:);
                            itest.Y1=Y1(testind,:);itest.Y2=Y2(testind,:);itest.Y3=Y3(testind,:);itest.Y4=Y4(testind,:);
                            agesTr=agesTrain(itrainind,:);snpsTr=SNPsTr(itrainind,:);
                            agesTe=agesTrain(itestind,:);snpsTest=SNPsTr(itestind,:);
                            %% reformat data
                            
                            dataTrain = ReformatData(snpsTr,agesTr,itrain.Y1,itrain.Y2,itrain.Y3,itrain.Y4);
                            iYtrain = dataTrain.Y; iAtrain = dataTrain.A; iXtrain = dataTrain.X; iXIndtrain=dataTrain.F;
                            dataTest = ReformatData(snpsTe,agesTe,itest.Y1,itest.Y2,itest.Y3,itest.Y4);
                            iYtest = dataTest.Y; iAtest = dataTest.A; iXtest = dataTest.X; iXIndtest=dataTest.F;
                            
                            % compute smoother matrix
                            
                            iStrain = GetSmoother(Atrain,Atrain,options);
                            iStest = GetSmoother(Atest,Atrain,options);
                            % get dimension
                            n=size(itrain.Y1,1);
                            p=size(Xtrain,2);
                            t=size(ages,2);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Run Algorithm %
                            
                            startTime = cputime;
                            
                            %% %choose lambda using a held out validation set
                            
                            tic
                            
                            % perform feature selection
                            tic
                            [snps,rawWeights,rawEstimates] = SelectGroups(itrain,iYtrain,iStrain,iXtrain,opts,p,options,iXIndtrain);
                            if (options.someVerbose)
                                fprintf('Selected %d SNPs\n',length(snps));
                            end
                            toc
                            % refit model with selected groups
                            tic
                            [refitWeights,refitEstimates] = EstimateFunctions(itrain,iYtrain,iStrain,iXtrain(:,snps),options,rawEstimates);
                            toc
                            %             %predict traits on validation set
                            tic
                            
                            iYpred = PredictEffect(itrain,iXtrain(:,snps),iStest,iXtest(:,snps),refitEstimates,rawEstimates);
                            toc
                            iV=rawEstimates.V;
                            iYtest1=test.Y1*V(:,1);
                            iYtest2=test.Y1*V(:,2);
                            iYtest3=test.Y1*V(:,3);
                            iYtest4=test.Y1*V(:,4);
                            data1= ReformatData(isnpsTest,ages,iYtest1,iYtest2,iYtest3,iYtest4);
                            iYtest=data1.Y;
                            res_kfold_RMSE(i)=sqrt(mean((iYtest-iYpred).^2));
                            
                        end
                        res_RMSE(j1,j2,j3,j4,j5)=mean(res_kfold_RMSE);
                    end
                end
            end
        end
    end
    ndim=size(res_RMSE);
    tempRMSE=10;
    for j1=1:ndim(1)
        for j2=1:ndim(2)
            for j3=1:ndim(3)
                for j4=1:ndim(4)
                    for j5=1:ndim(5)
                        if  RMSE1(j1,j2,j3,j4,j5)<tempRMSE
                            tempRMSE=res_RMSE(j1,j2,j3,j4,j5);
                            para=[j1,j2,j3,j4,j5];
                        end
                    end
                end
            end
        end
    end
    opts.lambdaVals.v2=lambdaVals1(para);
    opts.lambdaVals.v2=lambdaVals1(para);
    opts.lambdaVals.v2=lambdaVals1(para);
    opts.lambdaVals.v2=lambdaVals1(para);
    opts.lambdaVals.v2=lambdaVals1(para);
    dataTrain = ReformatData(snpsTr,agesTr,itrain.Y1,itrain.Y2,itrain.Y3,itrain.Y4);
    Ytrain = dataTrain.Y; Atrain = dataTrain.A; Xtrain = dataTrain.X; XIndtrain=dataTrain.F;
    dataTest = ReformatData(snpsTe,agesTe,itest.Y1,itest.Y2,itest.Y3,itest.Y4);
    Ytest = dataTest.Y; Atest = dataTest.A; Xtest = dataTest.X; XIndtest=dataTest.F;
    
    % compute smoother matrix
    
    Strain = GetSmoother(Atrain,Atrain,options);
    Stest = GetSmoother(Atest,Atrain,options);
    [snps,rawWeights,rawEstimates] = SelectGroups(train,Ytrain,Strain,Xtrain,opts,p,options,XIndtrain);
    if (options.someVerbose)
        fprintf('Selected %d SNPs\n',length(snps));
    end
    % refit model with selected groups
    [refitWeights,refitEstimates] = EstimateFunctions(train,Ytrain,Strain,Xtrain(:,snps),options,rawEstimates);
    
    %             %predict traits on validation set
    [Ypred,fhattest] = PredictEffect(train,Xtrain(:,snps),Stest,Xtest(:,snps),refitEstimates,rawEstimates);
    V=rawEstimates.V;
    Ytest1=test.Y1*V(:,1);
    Ytest2=test.Y1*V(:,2);
    Ytest3=test.Y1*V(:,3);
    Ytest4=test.Y1*V(:,4);
    data1= ReformatData(snpsTest,agesTest,Ytest1,Ytest2,Ytest3,Ytest4);
    iYtest=data1.Y;
    testRMSE(k)=sqrt(mean((iYtest-iYpred).^2));
end
RMSE=[mean(testRMSE) std(testRMSE)]
endTime = cputime;
runTime = endTime - startTime;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Results %
save(fullfile(outputDir,[options.outputFile(1:end-4) '.mat']),'gwasResults','runTime','RMSE');



