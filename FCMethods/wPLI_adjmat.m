function wPLI_adjmat=wPLI_adjmat(data)
[reg,~]=size(data);
wPLI_adjmat= zeros(reg,reg);
for i=1:reg
    k=1+i;
    for j=k:reg
        wPLI_adjmat(i,j)=wPLI_compute(data(i,:),data(j,:));        
        wPLI_adjmat(j,i)=wPLI_adjmat(i,j);
    end
end
end