function [G,S] = synthG(synsteps,Ht,Hc,Hp,X)
% [G,S] = synthG(synsteps,Ht*,Hc*,Hp*,X*)
% 
% synthG calculates the 1-dimensional synthesis gain for a unit excitation
% in a dwt band.
% The sequence of dwt synthesis and externally specified filters used in
% the synthesis steps must be specified through the synsteps argument,
% which is a string composed of 'L', 'H', 'T', 'C' and 'P' chars, asserted
% in analysis order.
% Multiple parallel synthesis sequences can be specified by comma
% separating within the control string. These parallel synthesis vectors
% will be summed.
% 
% The synthesis vector is computed by working backwards through the
% synsteps string:
% 'L': interleave with zeros then convolve with the lowpass 9x7 dwt
% 'H': interleave with zeros then convolve with the highpass 9x7 dwt
% 'T': convolve with -Ht
% 'C': convolve with Hc
% 'P': convolve with Hp
%
% If supplied, X is used as the excitation vector instead of the unit
% impulse. This mode can be used if synthG is called to continue expansion
% of an existing synthesis vector.

if any(synsteps=='X') % do nothing command
    S=0; G=0; return;
end
if all(ismember(synsteps,'LHTCP,'))
    N = length(synsteps);
else error('Unknown control character entered')
end
if (nargin<2)
    Ht=1; end
if (nargin<3)
    Hc = fliplr(Ht); end
if (nargin<4)
    Hp=Hc; end
if (nargin<5)
    X=1; end

% dwt 9/7 synthesis vectors
low = [-2*0.045636 -2*0.028772 2*0.295636 2*0.557543 2*0.295636 -2*0.028772 -2*0.045636];
high = [2*0.026749 2*0.016864 -2*0.078223 -2*0.266864 ...
    2*0.602949 -2*0.266864 -2*0.078223 2*0.016864 2*0.026749];

m=1; Y{m}=X;
for n = N:-1:1
    T = zeros(1,2*length(Y{m})-1);
    T(1:2:length(T)) = Y{m}; % interleaved with zeros
    switch synsteps(n)
        case 'L' % dwt; expand with zeros then convolve
            Y{m} = conv(T,low);
        case 'H'
            Y{m} = conv(T,high);
        case 'T' % packet lifting filters; just convolve
            Y{m} = -conv(Y{m},Ht); % subtracted during synthesis
        case 'C'
            Y{m} = conv(Y{m},Hc);
        case 'P'
            Y{m} = conv(Y{m},Hp);
        case ',' % start another synthesis path
            m=m+1; Y{m}=X;
    end
end

% combine synthesis paths
S=Y{1};
for n=2:m
    S = sumvectors(S,Y{n});
end
G = sum(S.^2);