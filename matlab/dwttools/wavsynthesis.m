function out = wavsynthesis(in,spacing,offset,g0,g1)
% out = wavsynthesis(in,spacing,offset,g0,g1)
%
% Performs the inverse of the transform wavanalysis

N = length(g1)-length(g0);
if N>0, g0 = [zeros(1,N/2) g0 zeros(1,N/2)];
else    g1 = [zeros(1,-N/2) g1 zeros(1,-N/2)];
end % zeropad shorter vector
N = max(length(g0),length(g1));

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

% Initialize array to hold output.
h = h + 2*overlap; % Treat extended input as input for now.
ovals = zeros(h+2*overlap,w);
% Apply low-pass synthesis basis vectors.
if mod(overlap,2)==0 % if overlap is even
    low_indices = 0:2:(h-1);
    high_indices = 1:2:(h-1);
else % symmetric extension with odd samples shifts indices
    low_indices = 1:2:(h-1);
    high_indices = 0:2:(h-1);
end
for n=1:N
    ovals(low_indices+n,:) = ovals(low_indices+n,:) + g0(n)*ivals(low_indices+1,:);
end
% Apply high-pass synthesis basis vectors.
for n=1:N
    ovals(high_indices+n,:) = ovals(high_indices+n,:) + g1(n)*ivals(high_indices+1,:);
end
% Overwrite relevant locations in `out' array with reconstructed samples
h = h - 2*overlap; % Restore length of active region.
out(indices,:) = ovals(1+2*overlap:1:(h+2*overlap),:);