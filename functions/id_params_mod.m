function id_params_list=id_params_mod(idn_mat,bands,sub)

    % id_params() compute the parameters of ID matrix namely Idiff, Iself and
    % Iothers.
    id_params_list=zeros(4,5);
    Iself=zeros(1,bands);
    Iothers=zeros(1,bands);
    Idiff=zeros(1,bands);
    srate=zeros(1,bands);
    mean_mat=zeros(sub,bands);
    mask_diag=logical(eye(sub));
    for k=1:bands
        for i=1:sub
            sum_mat=sum(idn_mat(i,:,k))+sum(idn_mat(:,i,k))-2*(idn_mat(i,i,k));
            mean_mat(i,k)=sum_mat/((sub*2)-1);
        end
        temp_mat=idn_mat(:,:,k);
        %id params
        Iself(1,k)=mean(temp_mat(mask_diag));
        Iothers(1,k)=mean(mean_mat(:,k));
        Idiff(1,k)=mean(temp_mat(mask_diag)-mean_mat(:,k));
        %success rate
        srate(1,k)=sucrate(idn_mat(:,:,k),sub);
        
        id_params_list(1,:)=Iself;
        id_params_list(2,:)=Iothers;
        id_params_list(3,:)=Idiff;
        id_params_list(4,:)=srate;        
    end
end