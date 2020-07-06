%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to create expanded data and collapsed data structs %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edata,cdata] = ReformatData(snps,ages,Y1,Y2,Y3,Y4)

% sizes
n = size(Y1,1);
p = size(snps,2);
t = 4;
q=size(Y1,2);

% structs
edata = [];
cdata = [];

% features
edata.F = zeros(n*t,3*p);
cdata.F = zeros(n,3*p);
for i = 1:n
    for k = 1:t
        ind = (i-1)*t + k;
        for g = 1:3
            edata.F(ind,g:3:end) = snps(i,:) == g-1;
            cdata.F(i,g:3:end) = snps(i,:) == g-1;
        end
    end
end

% snps
edata.X = zeros(n*t,p);
cdata.X = snps;
for i = 1:n
    for k = 1:t
        ind = (i-1)*t + k;
        edata.X(ind,:) = snps(i,:);
    end
end

% ages
edata.A = zeros(n*t,1);
cdata.A = mean(ages,2);
for i = 1:n
    for k = 1:t
        ind = (i-1)*t + k;
        edata.A(ind) = ages(i,k);
    end
end

% phenotypes
edata.Y = zeros(n*t,q);
for j=1:n
    ind=(j-1)*t+(1:t);
    edata.Y(ind,:)=[Y1(j,:);Y2(j,:);Y3(j,:);Y4(j,:)];
end

% identities
edata.I = zeros(n*t,1);
cdata.I = (1:n)';
for i = 1:n
    for k = 1:t
        ind = (i-1)*t + k;
        edata.I(ind) = i;
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%