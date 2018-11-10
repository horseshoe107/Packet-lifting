function out = analysis_5_3(in,spacing,offset)
%     [out] = analysis_5_3(in,spacing,offset)
% Performs 1D analysis using the 5/3 transform. Each column of the `in'
% array is processed independently.  Within each column vector, the samples
% which input data is considered to belong to a sub-sequence of the
% sequence represented by the entire column.  The first entry in the
% sub-sequence starts `offset' entries into the vector and the elements of
% the sub-sequence are separated by `spacing' entries within the vector.
% The returned `out' array is identical to `in', except that the relevant
% input samples are overwritten with the interleaved sequence of subband
% samples:
% low-pass samples in the 1'st, 3'rd, ... entries;
% high-pass samples in the 2'nd, 4'th, ... entries.
length = size(in,1); % number of rows
nvecs = size(in,2);
out = in;
indices = (1+offset):spacing:length;
ivals = [zeros(2,nvecs);in(indices,:);zeros(2,nvecs)];
length = size(ivals,1)-4;
% Perform symmetric extension.  Should work for inputs
% down to length 1.
for i=1:2
   ivals(3-i,:) = ivals(3+i,:);
   ivals(2+length+i,:) = ivals(2+length-i,:);
end
% Generate low-pass samples
low_indices = 0:2:(length-1);
low = 0.75*ivals(low_indices+3,:) ...
   + 0.25*ivals(low_indices+2,:) ...
   + 0.25*ivals(low_indices+4,:) ...
   - 0.125*ivals(low_indices+1,:) ...
   - 0.125*ivals(low_indices+5,:);
% Generate high-pass samples
high_indices = 1:2:(length-1);
high = 0.5*ivals(high_indices+3,:) ...
     - 0.25*ivals(high_indices+2,:) ...
     - 0.25*ivals(high_indices+4,:);
% Overwrite relevant locations in `out' array with
% interleaved subband samples.
interleaved = zeros(length,nvecs);
interleaved(1:2:length,:) = low;
interleaved(2:2:length,:) = high;
out(indices,:) = interleaved;