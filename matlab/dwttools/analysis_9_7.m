function out = analysis_9_7(in,spacing,offset)
% out = analysis_9_7(in,spacing,offset)
%
% Performs 1D analysis using the 9/7 transform. Each column of the `in'
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
length = size(in,1);
nvecs = size(in,2);
out = in;
indices = (1+offset):spacing:length;
ivals = [zeros(4,nvecs);in(indices,:);zeros(4,nvecs)];
length = size(ivals,1)-8;
% Perform symmetric extension.  Should work for inputs down to length 1.
for i=1:4
   ivals(5-i,:) = ivals(5+i,:);
   ivals(4+length+i,:) = ivals(4+length-i,:);
end
% Generate low-pass samples
low_indices = 0:2:(length-1);
low = 0.602949*ivals(low_indices+5,:) ...
    + 0.266864*ivals(low_indices+4,:) ...
    + 0.266864*ivals(low_indices+6,:) ...
    - 0.078223*ivals(low_indices+3,:) ...
    - 0.078223*ivals(low_indices+7,:) ...
    - 0.016864*ivals(low_indices+2,:) ...
    - 0.016864*ivals(low_indices+8,:) ...
    + 0.026749*ivals(low_indices+1,:) ...
    + 0.026749*ivals(low_indices+9,:);
% Generate high-pass samples
high_indices = 1:2:(length-1);
high = 0.557543*ivals(high_indices+5,:) ...
     - 0.295636*ivals(high_indices+4,:) ...
     - 0.295636*ivals(high_indices+6,:) ...
     - 0.028772*ivals(high_indices+3,:) ...
     - 0.028772*ivals(high_indices+7,:) ...
     + 0.045636*ivals(high_indices+2,:) ...
     + 0.045636*ivals(high_indices+8,:);
% Overwrite relevant locations in `out' array with interleaved subband samples.
interleaved = zeros(length,nvecs);
interleaved(1:2:length,:) = low;
interleaved(2:2:length,:) = high;
out(indices,:) = interleaved;