function plot_graph(ID_mat,bands,sub, method_n)
    %% [Computing ID parameters]
    % Computing ID matrices parameters: Idiff, Iself and Iothers 

    % id_params() computes Idiff, Iself, and Iothers
    [Idiff, ~, ~]=id_params(ID_mat,bands,sub);

    % sucrate() computes success rate
    sr=zeros(5,1);
    for j=1:bands
        sr(j,1)=sucrate(ID_mat(:,:,j),sub);
    end
    %% [Plotting]
    % Plots the ID matrices for the 2 frequency bands (alpha and beta) and displays the ID
    % parameters for each ID matrix.

    figure;
    subplot(2,1,1)
    lablx= 'Session 1';
    lably= 'Session 2';
    imagesc(ID_mat(:,:,3)),xlabel(lablx), ylabel(lably),...
        colorbar, title([method_n,' | Alpha band',' | Idiff=', num2str(Idiff(3)),' | SR=',num2str(sr(3)), '%']) 
    axis square;
    subplot(2,1,2)
    imagesc(ID_mat(:,:,4)),xlabel(lablx), ylabel(lably),...
        colorbar, title([method_n,' | Beta band',' | Idiff=', num2str(Idiff(4)),' | SR=',num2str(sr(4)), '%']) 
    colormap(jet);
    axis square;
end
