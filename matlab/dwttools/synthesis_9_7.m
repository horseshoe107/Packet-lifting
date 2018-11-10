function out = synthesis_9_7(in,spacing,offset)
% out = synthesis_9_7(in,spacing,offset)
% 
% Inverts the transform applied by function "analysis_9_7".
% This function has been rewritten in terms of predefined matrix
% operations. Profiler reports a speedup of ~x2 as a result.
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
% Initialise output array and define basis vectors
length = length + 8; % Treat extended input as input for now.
% NEW CODE BELOW
low_indices = 1:2:length;
high_indices = 2:2:length;
ovals = zeros(length,nvecs);
% Apply low pass synthesis
lowsynth = 2*[-0.045636 -0.028772 0.295636 0.557543 0.295636 -0.028772 -0.045636];
ivals2 = ivals; ivals2(high_indices,:)=0;
ovals = ovals+conv2(lowsynth,1,ivals2,'same');
% Apply high pass synthesis
highsynth = 2*[0.026749 0.016864 -0.078223 -0.266864...
    0.602949 -0.266864 -0.078223 0.016864 0.026749];
ivals2 = ivals; ivals2(low_indices,:)=0;
ovals = ovals+conv2(highsynth,1,ivals2,'same');
% Overwrite relevant locations in `out' array with reconstructed samples
length = length - 8; % Restore length of active region.
out(indices,:) = ovals(5:(length+4),:);

% DEPRECATED CODE
% ovals = zeros(length+8,nvecs); % Output has another 4 extra samples before and after active region.
% % Apply low-pass synthesis basis vectors.
% low_indices = 0:2:(length-1);
% ovals(low_indices+5,:) = ovals(low_indices+5,:) + 2*0.557543*ivals(low_indices+1,:);
% ovals(low_indices+4,:) = ovals(low_indices+4,:) + 2*0.295636*ivals(low_indices+1,:);
% ovals(low_indices+6,:) = ovals(low_indices+6,:) + 2*0.295636*ivals(low_indices+1,:);
% ovals(low_indices+3,:) = ovals(low_indices+3,:) - 2*0.028772*ivals(low_indices+1,:);
% ovals(low_indices+7,:) = ovals(low_indices+7,:) - 2*0.028772*ivals(low_indices+1,:);
% ovals(low_indices+2,:) = ovals(low_indices+2,:) - 2*0.045636*ivals(low_indices+1,:);
% ovals(low_indices+8,:) = ovals(low_indices+8,:) - 2*0.045636*ivals(low_indices+1,:);
% % Apply high-pass synthesis basis vectors.
% high_indices = 1:2:(length-1);
% ovals(high_indices+5,:) = ovals(high_indices+5,:) + 2*0.602949*ivals(high_indices+1,:);
% ovals(high_indices+4,:) = ovals(high_indices+4,:) - 2*0.266864*ivals(high_indices+1,:);
% ovals(high_indices+6,:) = ovals(high_indices+6,:) - 2*0.266864*ivals(high_indices+1,:);
% ovals(high_indices+3,:) = ovals(high_indices+3,:) - 2*0.078223*ivals(high_indices+1,:);
% ovals(high_indices+7,:) = ovals(high_indices+7,:) - 2*0.078223*ivals(high_indices+1,:);
% ovals(high_indices+2,:) = ovals(high_indices+2,:) + 2*0.016864*ivals(high_indices+1,:);
% ovals(high_indices+8,:) = ovals(high_indices+8,:) + 2*0.016864*ivals(high_indices+1,:);
% ovals(high_indices+1,:) = ovals(high_indices+1,:) + 2*0.026749*ivals(high_indices+1,:);
% ovals(high_indices+9,:) = ovals(high_indices+9,:) + 2*0.026749*ivals(high_indices+1,:);
% % Overwrite relevant locations in `out' array with reconstructed samples
% length = length - 8; % Restore length of active region.
% out(indices,:) = ovals(9:1:(length+8),:);