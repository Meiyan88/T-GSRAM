%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate f %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edata= reformatf(f)

% sizes
t = 4;
n = size(f,1)/t;
q=size(f,2);

% structs
edata = [];

% phenotypes
edata.f1 = zeros(n,q);edata.f2 = zeros(n,q);
edata.f3 = zeros(n,q);edata.f4 = zeros(n,q);
for j=1:n
    ind=(j-1)*t+(1:t);
    edata.f1(j,:)= f(ind(1),:);
    edata.f2(j,:)= f(ind(2),:);
    edata.f3(j,:)= f(ind(3),:);
    edata.f4(j,:)= f(ind(4),:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%