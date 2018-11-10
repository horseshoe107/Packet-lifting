function out = wavanalysis(in,spacing,offset,h0,h1)
% out = wavanalysis(in,spacing,offset,h0,h1)
%
% Analyse the columns of `in' with externally specified wavelet analysis
% kernels, h0 and h1. Both h0 and h1 are assumed to have support of odd
% length.

N = length(h1)-length(h0);
if N>0, h0 = [zeros(1,N/2) h0 zeros(1,N/2)];
else    h1 = [zeros(1,-N/2) h1 zeros(1,-N/2)];
end % zeropad shorter vector
N = max(length(h0),length(h1));

overlap = (N-1)/2; % N must be odd; undefined otherwise
[h,w] = size(in);
out = in;
indices = (1+offset):spacing:h;
ivals = [zeros(overlap,w);in(indices,:);zeros(overlap,w)];
h = size(ivals,1)-2*overlap;
% Perform symmetric extension.  Should work for inputs down to length 1.
for i=1:overlap
   ivals(overlap+1-i,:) = ivals(overlap+1+i,:);
   ivals(overlap+h+i,:) = ivals(overlap+h-i,:);
end

% Generate low-pass samples
low_indices = 0:2:(h-1);
low = zeros(ceil(h/2),w);
for n=1:N
    low = low + h0(n)*ivals(low_indices+n,:);
end
% Generate high-pass samples
high_indices = 1:2:(h-1);
high = zeros(floor(h/2),w);
for n=1:N
    high = high + h1(n)*ivals(high_indices+n,:);
end
% Overwrite relevant locations in `out' array with interleaved subband samples.
interleaved = zeros(h,w);
interleaved(1:2:h,:) = low;
interleaved(2:2:h,:) = high;
out(indices,:) = interleaved;