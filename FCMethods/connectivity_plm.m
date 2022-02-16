function [p] = connectivity_plm(input, B, fs)

% CONNECTIVITY_PLM computes the phase linearity measurement from a cell
% array of time-domain data, where each cell is an epoch. The methodology
% is described in Baselice et al. "Phase Linearity Measurement: a novel index
% for brain functional connectivity", IEEE TMI, 2018.
%
% Use as
%   [p] = ft_connectivity_plm(input, ...)
%
% The input data is a matrix of of nchan x nsamples signals (one epoch).
%
% Additional input argument
%   B	=	scalar, half-bandwidth parameter: the frequency range
%			across which to integrate (default=1 Hz)
%   fs  =       sampling frequency, needed to convert bandwidth to number of bins
%
% The output p contains the phase linearity measurement in the [0, 1] range.
% The output p is organized as a matrix of nchan x  nchan values.

input = hilbert(input.').';

nchan=size(input,1);
trial_length=size(input,2);
ph_min=0.1;        % Eps of Eq.(17) of the manuscript
f=(fs/trial_length)*(0:(trial_length-1));
f_integr=(abs(f)<B) | (abs(f-fs)<B);
p=zeros(nchan, nchan); % output variable

% Loop over channel pairs
for kchan1=1:(nchan-1)
    for kchan2=(kchan1+1):nchan
        
        % FFT of normalized interferometrix signal
        temp = fft( exp(1i*angle(input(kchan1,:))) .* exp(-1i*angle(input(kchan2,:))) );
        
        % FFT in zero ~ check phase of FFT in 0
        temp(1) = temp(1) .* (abs(angle(temp(1))) > ph_min);  % Volume conduction suppression
        
        % Take the square modulus of the FFT (energy spectral density)
        temp=(abs(temp)).^2;
        
        % Normalize the integral of the energy spectral density in the predefined band over the whole energy spectral density integral 
        p_temp=sum(temp(f_integr))./sum(temp);
        
        % PLM is symmetric
        p(kchan1, kchan2)=p_temp;
        p(kchan2, kchan1)=p_temp;
    end
end
