function out = synthesis_5_3(in,spacing,offset)
%   [out] = synthesis_5_3(in,spacing,offset)
% Inverts the transform applied by function
% "analysis_5_3".
length = size(in,1);
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
% Initialize array to hold output.
length = length + 4; % Treat extended input as input for now.
ovals = zeros(length+4,nvecs); % Output has another 2 extra samples before and after active region.
% Apply low-pass synthesis basis vectors.
low_indices = 0:2:(length-1);
ovals(low_indices+3,:) = ovals(low_indices+3,:) + 2*0.5*ivals(low_indices+1,:);
ovals(low_indices+2,:) = ovals(low_indices+2,:) + 2*0.25*ivals(low_indices+1,:);
ovals(low_indices+4,:) = ovals(low_indices+4,:) + 2*0.25*ivals(low_indices+1,:);
% Apply high-pass synthesis basis vectors.
high_indices = 1:2:(length-1);
ovals(high_indices+3,:) = ovals(high_indices+3,:) + 2*0.75*ivals(high_indices+1,:);
ovals(high_indices+2,:) = ovals(high_indices+2,:) - 2*0.25*ivals(high_indices+1,:);
ovals(high_indices+4,:) = ovals(high_indices+4,:) - 2*0.25*ivals(high_indices+1,:);
ovals(high_indices+1,:) = ovals(high_indices+1,:) - 2*0.125*ivals(high_indices+1,:);
ovals(high_indices+5,:) = ovals(high_indices+5,:) - 2*0.125*ivals(high_indices+1,:);
% Overwrite relevant locations in `out' array with
% reconstructed samples
length = length - 4; % Restore length of active region.
out(indices,:) = ovals(5:1:(length+4),:);