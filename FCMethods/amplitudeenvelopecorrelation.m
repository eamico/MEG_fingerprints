function [AEC, AECc] = amplitudeenvelopecorrelation(data)
% Compute leakage-corrected and uncorrected Amplitude Envelop Correlation (AECs)
% Arjan Hillebrand
% VUmc, Amsterdam

% Number of channels / ROIs
M = size(data,2);
% Number of time points
T = size(data,1);


% Initilaize output matrix
AECc = zeros(M,M);
AEC = zeros(M,M);

% Analytic signal input time series 
% If data is a matrix, then hilbert finds the analytic signal corresponding to each column.
htx = hilbert(data);

% Cut initial and end points from hilbert transform
cut = 10; 
htx(end-cut:end,:) = [];
htx(1:cut,:) = [];

% Signal magnitude (or envelope)
envelope_x = abs(htx);
AEC = corr(envelope_x);

% Loop over channels / ROIs [i]
for i = 1:M

    % Time series of channel/ROI i (x)
    x = data(:,i);

    % Loop over channels / ROIs [j]
    for j = 1:M

        if i~=j

            % Time series of channel/ROI j (y)
            y = data(:,j);
            
            % Orthogonolize signal y with respcted to signal x
            % (pairwise leakage correction)
            [~,~,r] = regress(y,x); 

            % Analytic signal of r
            hty = hilbert(r);

            % Discard first and last time points of hilbert-transformed signal    
            hty(end-cut:end) = [];
            hty(1:cut) = [];

            % Signal magnitude (or envelope)
            envelope_y = abs(hty);   

            % Correlation between envelopes
            AECc(i,j) = corr(envelope_x(:,i),envelope_y);

        end

    end

end

% Average between the 'two directions' leakage correction
AECc = (AECc + AECc') ./ 2;

