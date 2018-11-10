function X = packliftanalysis(X,Ht,Hc,Hp,mode,Tgrid)
% function X = packliftanalysis(X,Ht,Hc,Hp,mode*,Tgrid*)
%                             * -  optional arguments
% applies antialiasing to the image X by filtering and lifting between
% packet subbands. by default, the transform will analyse the baseband
% image down to a packet structure, perform packet lifting, then synthesise
% back up to the usual, interleaved mallat (LL,HL,LH,HH) decomposition.
% 
% Ht, Hc and Hp are the tranfer, cancellation, and precompensate/update
% filters respectively. mode is a control string that may contain any of
% the following command characters:
% 'S' - assumes the data is already analysed down to the required packet
%       structure and skips all ordinary dwt steps, including synthesis.
%       don't set this if X is an unmodified baseband image.
% 'V' - performs packet lifting in the vertical direction only.
% 'P' - Hp will be interpreted as a precompensation filter.
% 'U' - Hp will be interpreted as an update filter. Note 'P' and 'U' cannot
%       both be set. If neither are set, Hp will be ignored and the
%       function will only perform two-step packet lifting.
%
% The optional argument Tgrid controls spatial adaptivity of the packet
% lifting. Tgrid should have dimensions h/4 by w/4 and each element should
% take a value in [0,1], where 0 indicates the packet lifting has been
% switched off at that location.
[h,w] = size(X);
if nargin<5,  mode = ''; end
if nargin<6,  Tgrid = ones(h/4,w/4); end
skipdwt = any(mode=='S');
verticalonly = any(mode=='V');
use5x3 = any(mode=='5');

if skipdwt==0
    if use5x3
        X = analysis_5_3(analysis_5_3(X',1,0)',1,0);
        X = analysis_5_3(X,2,0);
        X = analysis_5_3(X,2,1);
        X = analysis_5_3(X',2,0)';
        X = analysis_5_3(X',2,1)';
    else
        X = analysis_9_7(analysis_9_7(X',1,0)',1,0);
        X = analysis_9_7(X,2,0);
        X = analysis_9_7(X,2,1);
        X = analysis_9_7(X',2,0)';
        X = analysis_9_7(X',2,1)';
    end
end

% define the swap bands
LL_HL = X(1:4:h,3:4:w);
LL_LH = X(3:4:h,1:4:w);
LL_HH = X(3:4:h,3:4:w);
HL_LL = X(1:4:h,2:4:w);
HL_LH = X(3:4:h,2:4:w);
LH_LL = X(2:4:h,1:4:w);
LH_HL = X(2:4:h,3:4:w);

% vertical first
[LL_LH,LH_LL] = packswap(LL_LH,LH_LL,Ht,Tgrid,[mode 'A'],fliplr(Hc),fliplr(Hp));
[LL_HH,LH_HL] = packswap(LL_HH,LH_HL,Ht,Tgrid,[mode 'A'],fliplr(Hc),fliplr(Hp));
% horizontal
if (verticalonly==0)
    [LL_HL,HL_LL] = packswap(LL_HL,HL_LL,Ht,Tgrid',[mode 'AH'],fliplr(Hc),fliplr(Hp));
    [LL_HH,HL_LH] = packswap(LL_HH,HL_LH,Ht,Tgrid',[mode 'AH'],fliplr(Hc),fliplr(Hp));
end

% copyback bands
X(1:4:h,3:4:w) = LL_HL;
X(3:4:h,1:4:w) = LL_LH;
X(3:4:h,3:4:w) = LL_HH;
X(1:4:h,2:4:w) = HL_LL;
X(3:4:h,2:4:w) = HL_LH;
X(2:4:h,1:4:w) = LH_LL;
X(2:4:h,3:4:w) = LH_HL;

if skipdwt==0
    if use5x3
        X = synthesis_5_3(X',2,1)';
        X = synthesis_5_3(X',2,0)';
        X = synthesis_5_3(X,2,1);
        X = synthesis_5_3(X,2,0);
    else
        X = synthesis_9_7(X',2,1)';
        X = synthesis_9_7(X',2,0)';
        X = synthesis_9_7(X,2,1);
        X = synthesis_9_7(X,2,0);
    end
end