% Permutation Testing (PT) framework to assess the statistical significance of 
% the differential identifiability and success rate scores. Script linked 
% with the methodology mentioned in the research article 
% Sareen et. al 2021, NeuroImage.
% 
% Requirements:
% ID_construct(), id_params(), sucrate(), and perm_con()

%% Loading data tensors
tensor1=tensor_sub_connmat_restsesh1;
tensor2=tensor_sub_connmat_restsesh2;

%% Generating true labels

% Defining number of subjects and bands
bands=5; sub=23;

% Generating ID matrix and parameters with true labels
[~,~,ID_mat]= ID_construct(tensor1,tensor2,bands,sub);

% Idiff
[Idiff_true,~,~]=id_params(ID_mat,bands,sub);
% SR
sr_true=zeros(5,1);
for j=1:bands
    sr_true(j,1)=sucrate(ID_mat(:,:,j),sub);
end
%% Permutations
perms=1000;
[perm_idiff,perm_sr]= perm_con(tensor1,tensor2,bands,sub,perms);

%% Computing statistical significance (p-values)  
%% M1: Conservative method
p_cons_idiff=zeros(5,1);
p_cons_sr=zeros(5,1);
for k=1:bands
    p_cons_idiff(k,1)=(length(find(abs(perm_idiff(:,k)) > abs(Idiff_true(k))))+1)/(perms+1);
    p_cons_sr(k,1)=(length(find(abs(perm_sr(:,k)) > abs(sr_true(k))))+1)/(perms+1);
end
%% M2: Approximate method
% Uses the p_cons_idiff and p_cons_sr from M1 method. 
p_appxm=zeros(5,1);
for k=1:5
    pp=p_cons_idiff(k);
    mt=nchoosek(46,23);
    b=sum(abs(perm_idiff(:,k))>=abs(Idiff_true(k)));
    gg=@(pp)binocdf(b,perms,pp);
    p_appxm(k,1)=pp-integral(gg,0,0.5/mt+1);
end
%% Plotting null distribution and significance (p-values)
method='AECc';
bandn={'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
band=3; % Alpha band
con='Rest Sesh1-Rest Sesh2';

% Idiff
figure; 
histogram(perm_idiff(:,band), 20);
hold on;
xlabel('Random differences');
ylabel('Count')
title(strcat('PT[IDIFF]', {' '}, 'M:',method, {' '}, 'B:', bandn{band}, {' '}, 'C:', con));
% To use M2 for p-value computation, update the variable name of the p value matrix 
od = plot(Idiff_true(band), 0, '*r', 'DisplayName', sprintf('Observed difference \np = %f', p_cons_idiff(band)));
legend(od);

% SR
figure;
histogram(perm_sr(:,band));
hold on;
xlabel('Random differences');
ylabel('Count')
title(strcat('PT[SR]', {' '}, 'M:',method, {' '}, 'B:', bandn{band}, {' '}, 'C:', con));
% To use M2 for p-value computation, update the variable name of the p value matrix 
od = plot(sr_true(band), 0, '*r', 'DisplayName', sprintf('Observed difference \np = %f', p_cons_sr(band)));
legend(od);