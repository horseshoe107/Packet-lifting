function Hu = packlift_preupdate(Ht,Hc,varargin)
% Hu = packlift_preupdate(Ht,Hc,gamma*)
%
% For a given Ht and Hc to be used in packet lifting, this function designs
% a pre-update filter Hu which minimises the energy resulting from
% synthesis of the quantisation noise in the donor and receiver packet
% subbands. The filter is found by minimising the objective function:
%     J = J1 + J2
% where J1 is the contribution from the donor region Q noise (assumed to be
% uniformly distributed with power density gamma1) and J2 is the
% contribution from the receiver region (Q noise gamma2).
%          /                                    /
%     J1 = | gamma1 |a(w)|^2 dw     and    J2 = | gamma2 |b(w)|^2 dw
%          /                                    /
% where a(w) = 1 - Hu(w)Ht(w),
%       b(w) = Hc(w) - Hu(w)res(w),
%       res(w) = Ht(w)Hc(w) - 1
% Refer to (bk2,pp136-9) for full derivation.
% 
% As Ht and Hc are supplied, and the target and weighting functions are
% based on combinations of these two filters, this function processes
% everything in the time domain.

N = length(Hc); % desired length of update filter
if nargin==3 % relative energy of quantisation noise for the band
    gamma1 = varargin{1};
else gamma1 = 1;
end
gamma2 = 1;
res = conv(Ht,Hc); res(ceil(length(res)/2)) = res(ceil(length(res)/2))-1;
% target function
d1 = fliplr(Ht);
d2 = conv(fliplr(res),Hc);
d = sumvectors(gamma1*d1,gamma2*d2)'; % col vector
% weighting function
r1 = conv(fliplr(Ht),Ht);
r2 = conv(fliplr(res),res);
r = sumvectors(gamma1*r1,gamma2*r2);
% form matrix equation and soln
d_overlap = length(d)-N;
d_leftskip = ceil(d_overlap/2);
d = d(d_leftskip+1:d_leftskip+N);
R = zeros(N);
r_centre = ceil(length(r)/2);
for y=1:N
    for x=1:N
        R(y,x) = r(r_centre+x-y);
    end
end
Hu = (inv(R)*d)';