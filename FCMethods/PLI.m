function PLI_mat= PLI(clean_data)
% Code adapted from the formulation explained in the original research 
% article titled: "Phase Lag Index: Assessment of Functional Connectivity 
% From Multi Channel EEG and MEG With Diminished Bias From Common Sources", 
% Cornelis J. Stam, Guido Nolte, and Andreas Daffertshofer, 
% Hum. Brain Mapp., 28 (2007), pp. 1178-1193.

    N = size(clean_data,1);
    PLI_mat=zeros(N,N);
    for i= 1:N
        for j= 1:N
            if i<j
                Channel_A=clean_data(i,:);
                Channel_B=clean_data(j,:);
                Hilbert_A=hilbert(Channel_A);
                Hilbert_B=hilbert(Channel_B);

                PhaseAng_A=angle(Hilbert_A);
                PhaseAng_B=angle(Hilbert_B);
                diff=wrapToPi(PhaseAng_A-PhaseAng_B);

                Size1=size(Channel_B);

                PLI_mat(i,j)=abs(sum(sign(diff)))/Size1(2);
            end               
        end
    end
   PLI_mat=PLI_mat+PLI_mat';
end