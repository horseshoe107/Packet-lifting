function X = fm2analysis(X,Ht,mode,varargin)
% function X = fm2analysis(X,Ht,mode*,Hc*,Hu/Hp*)
%                             * -  optional arguments
%
% fm2analysis is similar to fm_analysis, but a further iteration of dwt 9/7
% analysis is performed. Content is swapped between the L1 and H0 packet
% bands. As a result, the effective passband of the LPF used to transport
% content is reduced by a factor of 2. The 9/7 analysis filters are used to
% create the subbands.
%
% Ht should be a low pass filter. Ht is used to isolate the high frequency
% (aliased) data in the L1 band for transfer. If Ht is an FIR filter, its
% coefficients must be entered as a row vector. If Ht is IIR, the recursion
% coefficients must be entered in the 2nd row.
%
% If Hc is not explicitly specified, its coefficients are identical to Ht.
% Hc filters H0 noncausally during the cancellation step. Note that the
% end-to-end filter, conv(Ht,fliplr(Hc)) must be zero phase. It must also
% have a passband gain of 1. If ~Hc is obtained by factorisation of a
% designed zero phase FIR filter (Ht * ~Hc), fm2analysis should be supplied
% with fliplr(~Hc).
%
% The 3rd optional argument is a mode select string. Setting it to ''
% results in default behaviour.
% 'V' - only analyses in the vertical direction
% 'S' - assumes Xin is already analysed into interleaved dwt samples 
% 'SS' - assumes Xin is already analysed to TWO levels, in both the L and H
%        bands. Interleaving is LL,HL,LH,HH,LL... Synthesis out of this
%        packet level is not performed
% 'U' - insert an update step (order T,C,U). Hu=Ht if not specified
% 'P' - insert a precompensate step (P,T,C). Hp=Hc if not specified
%
% Hu/Hp is used in the update or pre-compensate step, depending on the
% mode. It may be explicitly specified as with Hc. Otherwise, it defaults
% to Ht (for 'U') or Hc (for 'P'). Note that for Hp, the same non-causal
% filtering concerns apply as for Hc.

[h,~]=size(X);

if nargin<3, mode=''; end % default empty string
horflag = all(mode~='V'); % set to 0 if 'V' is asserted
skipdwt = sum(mode=='S');
update  = any(mode=='U');
precomp = any(mode=='P');
if (update && precomp)
    error('3rd lifting step cannot be both a preupdate and an update!')
end

HtN=Ht(1,:); % extract transfer filter coefficients
if (size(Ht,1)>=2), HtD=Ht(2,:); % IIR
else HtD=1; end % FIR
% padding introduced by transfer filter
lpad=ceil(length(HtN)/2);
rpad=floor(length(HtN)/2);
if nargin>=4 % direct specification of the reverse filter
    HcN = varargin{1}(1,:);
    if (size(varargin{1},1)>=2), HcD=varargin{1}(2,:); % IIR
    else HcD=1; end
else HcN=HtN; HcD=HtD;
end
if nargin>=5 % direct specification of update filter
    HuN = varargin{2}(1,:);
    if (size(varargin{2},1)>=2), HuD=varargin{2}(2,:);
    else HuD=1; end
else
    if update
        HuN=HtN; HuD=HtD;
    else % assume precompensate
        HuN=HcN; HuD=HcD;
    end
end
% boundary extension: 0 for zero extension, -1 for zero incursion
bext = -1; % temporarily hardcoded

for i=0:horflag % iterate twice by default
    if horflag, X = X'; end % switch between hor/vert analyses
    
    % define the swap bands
    if (skipdwt==0), X=analysis_9_7(X,1,0); end % dwt analysis; skip for 'S'
    L = X(1:2:h,:);
    H = X(2:2:h,:);
    if (skipdwt<2) % further packet analysis; skip for 'SS'
        L = analysis_9_7(L,1,0);
        H = analysis_9_7(H,1,0);
    end % 'SS' flag: assume analysis already done externally
    Lswap = L(2:2:size(L,1),:); % LH
    Hswap = H(1:2:size(H,1),:); % HL
    
    % when H0~=L1 in size, we must conduct appropriate boundary interpolation
    boundH0=[1 , size(Hswap,1)]; % extend on RHS  (see bk2,p131)
    boundL1=[size(Hswap,1)-size(Lswap,1)+1 , size(Hswap,1)]; % restrict LHS noncausally
    if bext==-1
        boundL1 = [boundL1 rpad lpad];
    end

    % lifting steps between swap packets
    if precomp
        Y = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL1));
        Lswap = Lswap - Y;
    end
    Y = filtextend(Lswap,HtN,HtD,bext,boundH0); % transfer step
    Hswap = Hswap + Y;
    Y = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL1)); % cancellation - noncausal filtering
    Lswap = Lswap - Y;
    if update
        Y = filtextend(Lswap,HuN,HuD,bext,boundH0);
        Hswap = Hswap + Y;
    end

    L(2:2:size(L,1),:)=Lswap;
    H(1:2:size(H,1),:)=Hswap;

    if (skipdwt<2) % resynthesise after swapping
        X(1:2:h,:) = synthesis_9_7(L,1,0);
        X(2:2:h,:) = synthesis_9_7(H,1,0);
    else % 'SS' flag: leave as packets
        X(1:2:h,:) = L;
        X(2:2:h,:) = H;
    end
end