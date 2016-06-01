
% code for creating adjacency matrix from BOLD timeseries, with edge
% weights defined by magnitude-squared wavelet cross coherence, for one
% time window. 
% "window_ts" is an n x t matrix containing the time series for each node
% in one time window. (n = # nodes, t = # time samples in window)
% uses Grinsted et al. "wtc" package (http://noc.ac.uk/using-science/crosswavelet-wavelet-coherence)

[n,t] = size(window_ts);
bandpass = [0.06,0.125]; % in Hz
TR = 2;                  % TR in seconds

adj = zeros(n,n);        % initialise adjacency matrix
for i=1:n
    adj(i,i) = 0;        % set diagonal entries to 0
    for j=i+1:n
        % find cross coherence for one pair of regions i,j:
        [Rsq,period] = wtc(window_ts(i,:)',window_ts(j,:)','mcc',0);
        freq = 1./(period*TR);
        % average over time and over frequencies of interest:
        adj(i,j) = mean(mean(Rsq((freq<bandpass(2)&freq>bandpass(1)),:)));
    end
end
adj = adj + adj'; % force symmetric matrix