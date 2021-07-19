function [Idiff, Iself, Iothers]=id_params(idn_mat,bands,sub)

    % id_params() compute the parameters of ID matrix namely I diff, Iself and
    % Iothers.
    Iself=zeros(bands,1);
    Iothers=zeros(bands,1);
    Idiff=zeros(bands,1);
    mean_mat=zeros(sub,bands);
    mask_diag=logical(eye(sub));
    for k=1:bands
        for i=1:sub
            sum_mat=sum(idn_mat(i,:,k))+sum(idn_mat(:,i,k))-2*(idn_mat(i,i,k));
            mean_mat(i,k)=sum_mat/((sub*2)-1);
        end
        temp_mat=idn_mat(:,:,k);
        Iself(k)=mean(temp_mat(mask_diag));
        Iothers(k)=mean(mean_mat(:,k));
        Idiff(k)=mean(temp_mat(mask_diag)-mean_mat(:,k));
    end
end