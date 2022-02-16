function PLV_mat=PLV(clean_data)
% Code adapted from the formulation explained in the original research 
% article titled: Measuring phase synchrony in brain signals, J.P. Lachaux, 
% E. Rodriguez, J. Martinerie, F.J. Varela, Hum. Brain Mapp., 8 (1999), 
% pp. 194-208.
    
    N= size(clean_data,1);
    PLV_mat=zeros(N,N);
    for i= 1:N
        for j= 1:N
            if i<j
                Channel_A=clean_data(i,:);
                Channel_B=clean_data(j,:);
                Hilbert_A=hilbert(Channel_A);
                Hilbert_B=hilbert(Channel_B);
    
                PhaseAng_A=angle(Hilbert_A); 
                PhaseAng_B=angle(Hilbert_B);
                Size=size(Channel_B);
    
                PLV_mat(j,i)=abs(sum(exp(1j*(PhaseAng_A-PhaseAng_B))))/Size(2);
            end
        end
    end
    PLV_mat = PLV_mat + PLV_mat';
end