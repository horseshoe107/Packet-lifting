function [Lswap,Hswap] = packswap(Lswap,Hswap,Ht,Tgrid,mode,Hc,Hu)
% function [Lswap,Hswap] = packswap(Lswap,Hswap,Ht,Tgrid,mode*,Hc*,Hu/Hp*)
%
% Performs a series of antialiasing lifting steps between the packet
% subbands Lswap and Hswap. Adaptive control of the antialiasing can be
% applied through the Tgrid argument. Tgrid must be of dimension (h/4,w/4)
% and its elements should take values in the range [0,1]. Each element
% weights the amount added to Hswap at the corresponding pixel location
% during the transfer step. Thus, Tgrid=1 allows usual antialiasing
% operation, while Tgrid=0 disables it at that location.
if nargin<5, mode='A'; end % default settings
analysis = any(mode=='A');
horizontal = any(mode=='H');
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
if nargin<6
    HcN=HtN; HcD=HtD;
else
    HcN=Hc(1,:);
    if size(Hc,1)>=2, HcD=Hc(2,:);
    else              HcD=1; end
end
if nargin<7
    if update, HuN=HtN; HuD=HtD;
    else HuN=HcN; HuD=HcD; end % assume precompensate
else % direct specification of update filter
    HuN = Hu(1,:);
    if (size(Hu,1)>=2), HuD=Hu(2,:);
    else HuD=1; end
end

% boundary extension: 0 for zero extension, -1 for zero incursion
bext = -1; % temporarily hardcoded

if horizontal
    Lswap = Lswap';
    Hswap = Hswap';
end

% when H0~=L1 in size, we must conduct appropriate boundary interpolation
boundH0=[1 , size(Hswap,1)]; % extend on RHS  (see bk2,p131)
boundL1=[size(Hswap,1)-size(Lswap,1)+1 , size(Hswap,1)]; % restrict LHS noncausally
if bext==-1
    boundL1 = [boundL1 rpad lpad];
end

% lifting steps between swap packets
if analysis==1
    if precomp
        Y = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL1));
        Lswap = Lswap - Y;
    end
    Y = filtextend(Lswap,HtN,HtD,bext,boundH0); % transfer step
    Hswap = Hswap + Y.*Tgrid;
    % cancellation - noncausal filtering
    Y = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL1));
    Lswap = Lswap - Y;
    if update
        Y = filtextend(Lswap,HuN,HuD,bext,boundH0);
        Hswap = Hswap + Y;
    end
else
    if update
        Y = filtextend(Lswap,HuN,HuD,bext,boundH0);
        Hswap = Hswap - Y;
    end
    Y = flipud(filtextend(flipud(Hswap),HcN,HcD,bext,boundL1)); % reverse cancellation
    Lswap = Lswap + Y;
    Y = filtextend(Lswap,HtN,HtD,bext,boundH0); % reverse transfer step
    Hswap = Hswap - Y.*Tgrid;
    if precomp
        Y = flipud(filtextend(flipud(Hswap),HuN,HuD,bext,boundL1));
        Lswap = Lswap + Y;
    end
end

if horizontal
    Lswap = Lswap';
    Hswap = Hswap';
end