function icc_gen(tensor1,tensor2,method_n)
    %% [ICC matrix computing]

    % Note- Estimating ICC matrices only for frequency bands of interest i.e. 
    % alpha (band=3) and beta (band=4). For delta, theta and gamma band change
    % band variable to 1, 2 and 5 respectively. 

    % Alpha band
    ICC_mat_alpha=zeros(148,148);
    band=3;
    for i=1:148
        for j=1:148
            if i<j
                temp1=squeeze(tensor1(i,j,band,:));
                temp2=squeeze(tensor2(i,j,band,:));
                temp_comp=[temp1, temp2];
                [r, ~, ~, ~, ~, ~, ~] = ICC(temp_comp,'1-1', 0.05, 0);
                ICC_mat_alpha(i,j)=r;
            end                
        end
    end
    ICC_mat_alpha=ICC_mat_alpha+ICC_mat_alpha';

    % Beta band
    ICC_mat_beta=zeros(148,148);
    band=4;
    for i=1:148
        for j=1:148
            if i<j
                temp1=squeeze(tensor1(i,j,band,:));
                temp2=squeeze(tensor2(i,j,band,:));
                temp_comp=[temp1, temp2];
                [r, ~, ~, ~, ~, ~, ~] = ICC(temp_comp,'1-1', 0.05, 0);
                ICC_mat_beta(i,j)=r;
            end                
        end
    end
    ICC_mat_beta=ICC_mat_beta+ICC_mat_beta';
    %% ICC plotting
    figure;
    band_n1='Alpha';    % Band name 1
    band_n2='Beta';     % Band name 2
    % Note- Preload the Yeo parcellation. "aparc_a2009_yeo_RS7_MEG.mat" for
    % successful visualization.
    load('aparc_a2009_yeo_RS7_MEG.mat')
    
    %[Alpha band]
    % Weighted (fully-connected) ICC matrix
    subplot(2,1,1);
    imagesc(ICC_mat_alpha(yeoOrder,yeoOrder)); axis square; colormap(jet); colorbar
    title(strcat('ICC',{' '},'M:',method_n,{' '},'Band:', band_n1,{' '})) 

    %[Beta band]
    % Weighted (fully-connected) ICC matrix
    subplot(2,1,2);
    imagesc(ICC_mat_beta(yeoOrder,yeoOrder)); axis square; colormap(jet); colorbar
    title(strcat('ICC',{' '},'M:',method_n,{' '},'Band:', band_n2,{' '}))
end