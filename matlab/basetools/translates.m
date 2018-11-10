function [Y] = translates(A,B,N,offset)
% [Y] = translates(A,B,N*,offset*)
%
% Finds the dot product between vectors A and B, and all the N-translates
% between them. If N is neglected or set to 0, the function will calculate
% just a single dot product.
%
% A and B must be equal in length, or differ by an even number of samples.
% The default dot product is taken with A and B positioned such that their
% centres align. If the offset argument is used, B is shifted relative to A
% - to the right for offset>0, and left for offset<0
%
% Note: abs(offset) must be less than N.

% check A and B are vectors
if (ndims(A)>2)||(ndims(B)>2)||all(size(A)~=1)||all(size(B)~=1)
    error('Input arguments must be vectors.');
end
% convert to row vectors if not in the right form
if size(A,1)~=1, A=A'; end
if size(B,1)~=1, B=B'; end

olap = length(A) - length(B);
if mod(olap,2)~=0
    error('A and B must differ in length by an even number of samples');
end
if olap>0 % A longer, expand B
    B = [zeros(1,olap/2) B zeros(1,olap/2)];
else % overlap<=0
    A = [zeros(1,-olap/2) A zeros(1,-olap/2)];
end
L = length(A); % length of each vector (shorter one has been padded)

if (nargin<3)||(N==0)
    N=0;
    maxN=0;
else
    % maximum number of steps vectors can be translated relative to each other
    maxN = floor(L/N);
end
if nargin<4, offset=0; end

Y = zeros(1,2*maxN+1); % initialise
for n=-maxN:maxN
    shift = n*N+offset;
    if shift>0 % move B to the right
        D = [zeros(1,shift) B];
        C = [A zeros(1,shift)];
    else % negative shift; move B to the left
        D = [B zeros(1,-shift)];
        C = [zeros(1,-shift) A];
    end
    Y(maxN+n+1)= C*D'; % dot product
end