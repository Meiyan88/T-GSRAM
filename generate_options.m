clear
clc
%generate options struct
    options = []; 

% file names & data options

options.outputFile = 'gwas_results.mat';    % gwas results file

% parameter selection options

    options.numSnps = 200;                       % number of snps

    options.lambdaVals.v = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5];            % values of lambda
    
% algorithmic options

    options.kernelType = 'gauss';               % kernel type

    options.bandScaler = 1.5;                   % bandwidth scaler

    options.maxIter = 100;                      % max number of iterations
    options.convgTol = 1e-5;                    % tolerance for stopping criterion
    
% display & saving options

    options.someVerbose = true;                 % display printouts for each lambda
    
    save options.mat options
