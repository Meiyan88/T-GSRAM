%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to calculate smoother matrix between Atest and Atrain %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = GetSmoother(Atest,Atrain,options)

if ~exist('options','var') || ~isfield(options,'kernelType')
    ktype = 'gauss';
else
    ktype = options.kernelType;
end
if ~exist('options','var') || ~isfield(options,'bandScaler')
    scaler = 1.0; %ԭʼֵΪ1.0
else
    scaler = options.bandScaler;
end

[ntest,~] = size(Atest);  
[ntrain,~] = size(Atrain);

% plug-in bandwidth
bandwidth = scaler * std(Atrain) * ntrain^(-1/5);

% compute (n by n) smoother matrix
S = abs(repmat(Atest,1,ntrain) - repmat(Atrain',ntest,1))./bandwidth;

switch ktype
    case 'gauss'
        % Gaussian kernel
        S = gaussianKernel(S);
    case 'tgauss'
        % truncated Gaussian kernel
        S = sparse(truncGaussianKernel(S));
    case 'box'
        % Box kernel
        S = sparse(boxKernel(S));
    case 'triangle'
        % Triangular kernel
        S = sparse(triangularKernel(S));
    case 'epan'
        % Epanechnikov quadratic kernel
        S = sparse(epanechnikovKernel(S));
end

% normalize the rows of smoother matrix
S = normalizeRows(S);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization %

% normalize each row of A so that the row sum is 1
function A = normalizeRows(A)
    nCol = size(A,2);
    A = A .* repmat(1./sum(A, 2),1,nCol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels %

% Gaussian kernel
function K = gaussianKernel(x)
    K = normpdf(x,0,1);
end

% Truncated Gaussian kernel
function K = truncGaussianKernel(x)
    K = exp(-(x.^2)./2);
    K = (K >= exp(-0.5)).*K;
end

% Box kernel
function K = boxKernel(x)
    K = 0.5.*(double(abs(x) <= 1));
end

% Triangular kernel
function K = triangularKernel(x)
    K = max(0,1 - abs(x));
end

% Epanechnikov Quadratic kernel
function K = epanechnikovKernel(x)
    K = max(0,0.75.*(1 - x.^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%