function [triumask_sesh1, triumask_sesh2,ID_mat]= ID_construct(tensor_sesh1,tensor_sesh2,bands,sub)
    % ID_construct() computes the 2D matrix for each band and for each
    % session that contains the elements of vectorized upper triangle 
    % matrix ((N^2-N)/2 i.e. 10878 elements) of each subject for each band.
    % Using these matrices for each session Identifibility matrix (ID_mat) is computed.
    
    % Computing vectorized upper triangle matrix from each connectivity
    % matrix of every subject and every band
    triumask_sesh1=zeros(10878,sub,bands);
    triumask_sesh2=zeros(10878,sub,bands);
    
    ID_mat=zeros(sub,sub,bands);    
    mask=triu(true(148),1);
    for i=1:sub
        for j=1:bands
            temp_sesh1=tensor_sesh1(:,:,j,i);
            temp_sesh2=tensor_sesh2(:,:,j,i);
            triumask_sesh1(:,i,j)=temp_sesh1(mask);       
            triumask_sesh2(:,i,j)=temp_sesh2(mask);        
        end
    end
    % Computing ID matrix 
    for k=1:bands
        ID_mat(:,:,k) = corr(triumask_sesh1(:,:,k),triumask_sesh2(:,:,k));
    end
end