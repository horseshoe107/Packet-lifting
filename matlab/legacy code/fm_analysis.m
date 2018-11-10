function Y = fm_analysis(X,Ht,mode,varargin)
% function Y = fm_analysis(X,Ht,mode*,Hc*,Hu/Hp*)
%                             * -  optional arguments
%
% Performs one level of 2 dimensional dwt analysis, with the frequency
% modulation feature. The 9/7 analysis filters are used to create the
% subbands.
%
% Ht should be a low pass filter. Ht is used to isolate the high frequency
% (aliased) data in the L band for transfer. If Ht is an FIR filter, its
% coefficients must be entered as a row vector. If Ht is IIR, the recursion
% coefficients must be entered in the 2nd row.
%
% If Hc is not explicitly specified, its coefficients are identical to Ht.
% Hc filters H noncausally during the cancellation step. Note that the
% end-to-end filter, conv(Ht,fliplr(Hc)) must be zero phase. It must also
% have a passband gain of 1. If ~Hc is obtained by factorisation of a
% designed zero phase FIR filter (Ht * ~Hc), fm_analysis should be supplied
% with fliplr(~Hc).

% The 3rd optional argument is a mode select string. Setting it to ''
% results in default behaviour.
% 'V' - only analyse in the vertical direction
% 'S' - assumes X is already analysed into interleaved dwt samples 
% 'U' - insert an update step (order T,C,U). Hu=Ht if not specified
% 'P' - insert a precompensate step (P,T,C). Hp=Hc if not specified
%
% Hu/Hp is used in the update or pre-compensate step, depending on the
% mode. It may be explicitly specified as with Hc. Otherwise, it defaults
% to Ht (for 'U') or Hc (for 'P'). Note that for Hp, the same non-causal
% filtering concerns apply as for Hc.
if nargin<3, mode=''; end % default empty string
horflag = all(mode~='V'); % set to 0 if 'V' is asserted
skipdwt = sum(mode=='S');
update  = sum(mode=='U');
precomp = sum(mode=='P');
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
    if skipdwt==0, X=analysis_9_7(X,1,0); end
    % define swap bands
    Lswap = X(1:2:size(X,1),:);
    Hswap = X(2:2:size(X,1),:);
    % spectrum reversal kernel
    K = repmat((-1).^(1:ceil(size(X,1)/2))',1,size(X,2));
    % boundary interpolation
    boundL = [(size(Hswap,1)-size(Lswap,1))+1,size(Hswap,1)];
    boundH = [1,size(Hswap,1)];
    if bext==-1
        boundL1 = [boundL1 rpad lpad];
    end

    % frequency modulated lifting steps between swap bands
    if precomp
        A = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL));
        A = A .* K;
        Lswap = Lswap - A;
    end
    A = Lswap .* K;
    A = filtextend(A,HtN,HtD,bext,boundH);
    Hswap = Hswap + A;
    A = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL));
    A = A .* K;
    Lswap = Lswap - A;
    if update
        A = Lswap .* K;
        A = filtextend(A,HuN,HuD,bext,boundL);
        Hswap = Hswap + A;
    end

    X(1:2:size(X,1),:) = Lswap;
    X(2:2:size(X,1),:) = Hswap;
end
Y=X;