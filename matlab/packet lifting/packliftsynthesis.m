function X = packliftsynthesis(X,Ht,Hc,Hp,mode,Tgrid)
% function X = packliftsynthesis(X,Ht,Hc,Hp,mode*,Tgrid*)
%                             * -  optional arguments
% reverses the antialiasing packet lifting transform implemented by
% packliftanalysis. For a description of the arguments, see the analysis
% function.
[h,w] = size(X);
if nargin<5   mode = ''; end
if nargin<6   Tgrid = ones(h/4,w/4); end
skipdwt = any(mode=='S');
verticalonly = any(mode=='V');

if skipdwt==0
    X = analysis_9_7(X,2,0);
    X = analysis_9_7(X,2,1);
    X = analysis_9_7(X',2,0)';
    X = analysis_9_7(X',2,1)';
end

% define the swap bands
LL_HL = X(1:4:h,3:4:w);
LL_LH = X(3:4:h,1:4:w);
LL_HH = X(3:4:h,3:4:w);
HL_LL = X(1:4:h,2:4:w);
HL_LH = X(3:4:h,2:4:w);
LH_LL = X(2:4:h,1:4:w);
LH_HL = X(2:4:h,3:4:w);

% synthesis: horizontal first
if (verticalonly==0)
    [LL_HL,HL_LL] = packswap(LL_HL,HL_LL,Ht,Tgrid',[mode 'H'],fliplr(Hc),fliplr(Hp));
    [LL_HH,HL_LH] = packswap(LL_HH,HL_LH,Ht,Tgrid',[mode 'H'],fliplr(Hc),fliplr(Hp));
end
% vertical
[LL_LH,LH_LL] = packswap(LL_LH,LH_LL,Ht,Tgrid,[mode ''],fliplr(Hc),fliplr(Hp));
[LL_HH,LH_HL] = packswap(LL_HH,LH_HL,Ht,Tgrid,[mode ''],fliplr(Hc),fliplr(Hp));

% copyback bands
X(1:4:h,3:4:w) = LL_HL;
X(3:4:h,1:4:w) = LL_LH;
X(3:4:h,3:4:w) = LL_HH;
X(1:4:h,2:4:w) = HL_LL;
X(3:4:h,2:4:w) = HL_LH;
X(2:4:h,1:4:w) = LH_LL;
X(2:4:h,3:4:w) = LH_HL;

if skipdwt==0
    X = synthesis_9_7(X',2,1)';
    X = synthesis_9_7(X',2,0)';
    X = synthesis_9_7(X,2,1);
    X = synthesis_9_7(X,2,0);
end
X = synthesis_9_7(synthesis_9_7(X',1,0)',1,0);