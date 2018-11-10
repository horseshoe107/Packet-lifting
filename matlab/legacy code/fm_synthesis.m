function Y = fm_synthesis(X,Ht,mode,varargin)
% function Y = fm_synthesis(X,Ht,mode*,Hc*,Hu/Hp*)
%                             * -  optional arguments
%
% fm_synthesis is the synthesis complementary function of fm_analysis.
% If Xin was produced by fm_analysis, using identical control arguments
% will ensure perfect reconstruction of the original image. For a
% description of the arguments and their usage, refer to fm_analysis.
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
    
    % frequency modulated lifting steps in reverse (synthesis) order
    if update
        A = Lswap .* K;
        A = filtextend(A,HuN,HuD,bext,boundH);
        Hswap = Hswap - A;
    end
    A = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL));
    A = A .* K;
    Lswap = Lswap + A;
    A = Lswap .* K;
    A = filtextend(A,HtN,HtD,bext,boundH);
    Hswap = Hswap - A;
    if precomp
        A = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL));
        A = A .* K;
        Lswap = Lswap + A;
    end

    X(1:2:size(X,1),:) = Lswap;
    X(2:2:size(X,1),:) = Hswap;
    if skipdwt==0, X=synthesis_9_7(X,1,0); end
end
Y=X;