% wPLI computation between two series x and y
% This function is called within main function wPLI_adjmat.m
function wPLI_val=wPLI_compute(x,y)
    [~,len]=size(x);
    A=fft(x);
    B=fft(y);

    Channel_A_con_B=A(1:fix(len/2)).*conj(B(1:fix(len/2)));
    common=imag(Channel_A_con_B);
    num_1_wPLI=abs(mean(common));
    dem_1_wPLI=mean(abs(common));
    wPLI_val=num_1_wPLI/dem_1_wPLI;

end
