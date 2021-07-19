function sr_score=sucrate(ID_mat,sub)
    score_row=0;
    score_col=0;
    for i=1:sub
        temp_row=ID_mat(i,:);
        [~,idx1]=max(temp_row);
        if idx1==i
            score_row=score_row+1;
        end
    end
    for j=1:sub
        temp_col=ID_mat(:,j);
        [~,idx2]=max(temp_col);
        if idx2==j
            score_col=score_col+1;
        end
    end
    avg_score=(score_row+score_col)/2;
    sr_score=(avg_score/sub)*100;
end