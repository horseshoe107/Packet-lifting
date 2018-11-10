function [H0,H1] = liftexpand(varargin)
% [H0,H1] = liftexpand(L1,L2,L3,...)
% [G0,G1] = liftexpand('S',L1,L2,L3,...)
%
% Computes the low-pass and high-pass analysis and synthesis vectors from
% the individual filters in each step of a lifting structure.
%
% The arguments to liftexpand should be vectors containing the filter
% coefficients for each lifting step. The lifting structure used as a
% foundation has a summer at each predict/update node; subtraction must be
% produced by negating the corresponding filter's coefficients.
%
% The input analysed as 2 polyphase components: X0 and X1, which are the
% even and odd samples of the input respectively. The subscript
% n=[1:nargin] in either polyphase component refers to the progression
% along the structure as it passes through successive lifting steps.
%
% The X0 and X1 components are derived completely separately. The
% nomenclature YNXM denotes the support over the XMth polyphase component,
% required to produce output at the YNth polyphase component.
% The analysis function for YN is the complete support over X, so our final
% step should be a simple interleaving of YNX0 and YNX1.
%
% usage example for dwt9x7:
% l1 = repmat(-1.586134342,1,2);
% l2 = repmat(-0.052980118,1,2);
% l3 = repmat(0.882911075,1,2);
% l4 = repmat(0.443506852,1,2);
% K = 1.230174105; K0=1/K; K1=K/2;
% [Y0 Y1]=liftexpand(l1,l2,l3,l4); Y0=Y0*K0; Y1=Y1*K1;

synthflag=0;
ident = class(varargin{1});
if length(ident)==4 && all(ident=='char')
    if any(varargin{1}=='S') % reverse steps and filters now subtract
        N=nargin-1; synthflag=1;
        for n=1:N  L{n}=-varargin{nargin+1-n}; end
    else % assume analysis - copy remaining arguments across
        N=nargin-1;
        for n=1:N  L{n}=varargin{n+1}; end
    end
else
    L=varargin; N=nargin;
end

X0{1}=1; % excitations produced by even input samples
for n=1:N
    X0{n+1}=conv(X0{n},L{n});
    if (n>=2) X0{n+1}=sumvectors(X0{n+1},X0{n-1},0); end
end
Y=X0{N}; T=zeros(1,2*length(Y)-1); T(1:2:length(T))=Y;
if mod(N,2)==0, Y1X0=T; else Y0X0=T; end
Y=X0{N+1}; T=zeros(1,2*length(Y)-1); T(1:2:length(T))=Y;
if mod(N,2)==0 Y0X0=T; else Y1X0=T; end

X1{1}=1; % produced from odd input samples
for n=2:N
    X1{n}=conv(X1{n-1},L{n});
    if (n>=3) X1{n}=sumvectors(X1{n},X1{n-2},0); end
end
Y=X1{N-1}; T=zeros(1,2*length(Y)-1); T(1:2:length(T))=Y;
if mod(N,2)==0, Y1X1=T; else Y0X1=T; end
Y=X1{N}; T=zeros(1,2*length(Y)-1); T(1:2:length(T))=Y;
if mod(N,2)==0, Y0X1=T; else Y1X1=T; end

if synthflag==0 % analysis: interleave support over input
    H0=sumvectors(Y0X0,Y0X1,0);
    H1=sumvectors(Y1X0,Y1X1,0);
else % synthesis: interleave output from shared input
    if mod(N,2)==0 % first lifting step takes G1->G0 instead
        H0=sumvectors(Y0X1,Y1X1,0);
        H1=sumvectors(Y0X0,Y1X0,0);
    else
        H0=sumvectors(Y0X0,Y1X0,0);
        H1=sumvectors(Y0X1,Y1X1,0);
    end
end