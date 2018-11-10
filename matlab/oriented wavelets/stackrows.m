function M = stackrows(varargin)
% M = stackrows(R1,R2,....Rn)
%
% Form matrix M with row elements R1,R2..Rn such that M(i,:) = Ri
% The row vectors may be inequal in length but must differ from each other
% by an even number of elements. Shorter rows will be zero padded
% symmetrically to match the longest row.
for i=1:nargin
    R{i} = varargin{i};
    Rlen(i) = length(R{i});
end
L = max(Rlen);
M = zeros(nargin,L);
olap = L - Rlen;
if any(mod(olap,2)~=0)
    error('All rows must differ in length by only an even number of samples');
end
for i=1:nargin
    M(i,:) = [zeros(1,olap(i)/2) R{i} zeros(1,olap(i)/2)];
end