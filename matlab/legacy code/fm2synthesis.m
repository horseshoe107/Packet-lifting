function X = fm2synthesis(X,Ht,mode,varargin)
% function X = fm2synthesis(X,Ht,mode*,Hc*,Hu/Hp*)
%                              * -  optional arguments
%
% fm2synthesis is the synthesis complementary function of fm2analysis.
% If X was produced by fm2analysis, using identical control arguments
% will ensure perfect reconstruction of the original image. For a
% description of the arguments and their usage, refer to fm2analysis.

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
    
    L = X(1:2:h,:);
    H = X(2:2:h,:);
    if (skipdwt<2) % further packet analysis; skip for 'SS'
        L = analysis_9_7(X(1:2:h,:),1,0);
        H = analysis_9_7(X(2:2:h,:),1,0);
    end % 'SS' flag: assume analysis already done externally
    Lswap = L(2:2:size(L,1),:);
    Hswap = H(1:2:size(H,1),:);
    
    % when H0~=L1 in size, we must be careful about boundary interpolation
    boundH0=[1 , size(Hswap,1)]; % extend on RHS     (see bk2,p131)
    boundL1=[size(Hswap,1)-size(Lswap,1)+1 , size(Hswap,1)]; % restrict LHS noncausally
    if bext==-1
        boundL1 = [boundL1 rpad lpad];
    end
    
    % lifting in reverse order between swap packets
    if update
        Y = filtextend(Lswap,HuN,HuD,bext,boundH0);
        Hswap = Hswap - Y;
    end
    Y = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL1));
    Lswap = Lswap + Y;
    Y = filtextend(Lswap,HtN,HtD,bext,boundH0);
    Hswap = Hswap - Y;
    if precomp
        Y = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL1));
        Lswap = Lswap + Y;
    end
    
    L(2:2:size(L,1),:)=Lswap;
    H(1:2:size(H,1),:)=Hswap;

    if (skipdwt<2) % resynthesise after swapping
        X(1:2:h,:) = synthesis_9_7(L,1,0);
        X(2:2:h,:) = synthesis_9_7(H,1,0);
    else % 'SS' flag: leave in packets
        X(1:2:h,:) = L;
        X(2:2:h,:) = H;
    end
    if (skipdwt==0), A=synthesis_9_7(A,1,0); end % final dwt synthesis; skip for 'S'
end