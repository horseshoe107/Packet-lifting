function [Y,G] = sumvectors(A,B,varargin)
% [Y,G] = sumvectors(A,B,offset*)
%
% Adds 2 row vectors A and B together. When A and B are of equal length,
% this operation is trivial. If they differ in length, zero padding is
% applied appropriately to ensure the centres of the vectors are aligned.
% Note this operation requires the difference in length to be an even
% number of samples; when this is not satisfied, the output is
% unpredictable.
%
% B may be optionally shifted to the right relative to A by offset samples.
%
% The energy of the sum is returned in G
if nargin==3
    offset = varargin{1};
else offset = 0;
end

N = length(A)-length(B);
halfN=floor(N/2); ceilN=ceil(N/2);
if halfN ~= ceilN
    warning('WORK:sumvectors','Vectors do not have even difference in length.');
end
shift = max(min(offset,abs(halfN)),-abs(halfN)); % limit shift to [-N/2,N/2]
if N>0 % A longer, expand and shift B
    B = [zeros(1,ceilN+shift) B zeros(1,halfN-shift)];
else % B longer, expand A and shift in OPPOSITE direction
    A = [zeros(1,-ceilN-shift) A zeros(1,-halfN+shift)];
end

overshift = offset - shift;
if overshift>0 % move B to the right
    B = [zeros(1,overshift) B];
    A = [A zeros(1,overshift)];
else
    B = [B zeros(1,-overshift)];
    A = [zeros(1,-overshift) A];
end
Y = A + B;
G = sum(Y.^2);