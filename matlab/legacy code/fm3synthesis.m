function Xout = fm3synthesis(Xin,spacing,offset,M,Hb,Hl,beta,varargin)
% function Xout = fm3synthesis(Xin,spacing,offset,F,modes*,Frev*,r*)
%                             * -  optional arguments
%
% fm3synthesis is the synthesis complementary function of fm3analysis.
% If Xin was produced by fm3analysis, using identical control arguments
% will ensure perfect reconstruction of the original image. For a
% description of the arguments and their usage, refer to fm3analysis.
[h,w]=size(Xin);
hindices = (1+offset):spacing:h;
windices = (1+offset):spacing:w;
A = Xin(hindices,windices);
horflag=1; skipdwt=0;
if nargin>=8
    mode = varargin{1};
    if (any(mode=='V')), horflag=0; end % only vertical analysis
    if (any(mode=='S')), skipdwt=1; end
end

if skipdwt==0, A = synthesis_9_7(synthesis_9_7(A',1,0)',1,0); end
for i=0:horflag % iterate twice by default
    [h w] = size(A);
    olines=1:2:w/2; elines=2:2:w/2;
    
    A = analysis_9_7(analysis_9_7(A,1,0)',1,0)';

    % rearrange columns of H band before DWT synthesis
    HLolines = 1:4:w; HLelines = 3:4:w;
    A(2:2:h,HLelines) = A(2:2:h,HLelines) + A(2:2:h,HLolines);
    A(2:2:h,HLolines) = A(2:2:h,HLolines) - A(2:2:h,HLelines);
    A(2:2:h,HLelines) = A(2:2:h,HLelines) + A(2:2:h,HLolines);

    A = synthesis_9_7(A,1,0);

    L = A(:,1:2:w); % horizontal L band
    % Vertical analysis only. Swap using horizontal L bands only.
    % 3 lifting steps
    L(:,elines) = L(:,elines) - filtnflip(L(:,olines),M,Hb,Hl,beta);
    L(:,olines) = L(:,olines) + filtnflip(L(:,elines),M,Hb,Hl,beta);
    L(:,elines) = L(:,elines) - filtnflip(L(:,olines),M,Hb,Hl,beta);
    A(:,1:2:w) = L; % horizontal L band

    A = synthesis_9_7(A',1,0)';
    if horflag, A = A'; end % switch between vert/hor syntheses
end

Xout = Xin; Xout(hindices,windices)=A;