function [perm_idiff, perm_sr]= perm_con(tensor_sesh1,tensor_sesh2,bands,sub,perms)

% Function to generate the non-parametric null distribution of differential
% identifiability and success rate scores

    triumask_sesh1=zeros(10878,sub,bands);
    triumask_sesh2=zeros(10878,sub,bands);
    
    % Computing vectorized upper triangle matrix from each connectivity
    % matrix of every subject and every band
    mask=triu(true(148),1);
    for i=1:sub
        for j=1:bands
            temp_sesh1=tensor_sesh1(:,:,j,i);
            temp_sesh2=tensor_sesh2(:,:,j,i);
            triumask_sesh1(:,i,j)=temp_sesh1(mask);       
            triumask_sesh2(:,i,j)=temp_sesh2(mask); 
        end
    end
    
    % Random shuffling of the columns to break the session1 and session2
    % associations between the same subject and generate a random (null)
    % distribution
    for k=1:bands
        
        % Concatenating the #subjects in sesh1 and sesh2
        full_sub=[squeeze(triumask_sesh1(:,:,k)) squeeze(triumask_sesh2(:,:,k))];
        %[~,subss]=size(full_sub);

        for i=1:perms
            % Generating a random vector of subject size i.e. #subjects 
            % in session1
            g=randperm(sub);
            
            % Concatenating the random distribution of #subjects in sesh1 
            % while keeping the #subjects in sesh2 in the orginal order. 
            % This breaks the subjectwise sesh1 and sesh2 associations.
            h=[g 24:46];

            % Note- h [[#subjects in session 1] [#subjects in sessions 2]]
            % Here we have 23 subjects in session 1 and session 2; Hence
            % h=[1:23 24:46]
            
            % Ordering the columns on the new (random) order i.e. h
            rand_full_sub=full_sub(:,h);
            % Compute the ID matrix between the new (random) order 
            ID_mat = corr(rand_full_sub(:,1:23),rand_full_sub(:,24:46));
            % Computing the ID params
            [Idiff,~,~]=id_params(ID_mat,1,sub);
            
            % Storing the random ID parameters to generate the
            % distribution for statiscal test.
            perm_idiff(i,k)=Idiff;
            perm_sr(i,k)=sucrate(ID_mat,sub);
        end
    end
end