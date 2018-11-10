function Xout = fm3analysis(Xin,spacing,offset,M,Hb,Hl,beta,varargin)
% function Xout = fm3analysis(Xin,spacing,offset,M,Hb,Hl,modes*)
%                             * -  optional arguments
%
% fm3analysis also conducts modulated lifting, but we now apply this
% *before* dwt analysis (9/7). It is not a 1D function. For example, when
% vertically filtering the function only antialiases the L horizontal band,
% but transfers content to both LH and HH.
%
% Content may be modulated by an arbitrary amount; however, if we choose to
% shift it by pi/2, the resulting I and Q components are symmetrical about
% pi/2 in the frequency domain. We can then subsample by 2 BEFORE filtering
% with F (h_t). This special case is used in this function; this preserves
% some backwards compliance with fm_analysis, as the filter specifications
% are passband ? and transition band ?.
% NOTE if the special case is not met, we must subsample by 2 AFTER
% filtering, changing the filter specifications such that passband and
% transition band have width ?/2.
%
% If Hb or Hl are FIR filter, they must be entered as row vectors of
% coefficients. If they are IIR, the recursion coefficients must be entered
% in the 2nd row.
%
% The 5th optional argument is a string which enables different modes.
% Setting it to '' has no effect.
% If the char 'V' appears in this string, fm3analysis will only operate in
% the vertical direction.
% If the char 'S' appears once, fm3analysis will skip the initial dwt
% analysis and assume the subset A is already in the form of L and H
% interleaved bands.
[h,w]=size(Xin);
hindices = (1+offset):spacing:h;
windices = (1+offset):spacing:w;
A = Xin(hindices,windices);
horflag=1; skipdwt=0;
if nargin>=8
    mode = varargin{1};
    if (any(mode=='V')), horflag=0; end % block horizontal analysis
    if (any(mode=='S')), skipdwt=1; end
end

for i=0:horflag % iterate twice by default
    if horflag, A = A'; end % switch between hor/vert analyses

    [h w] = size(A);
    olines=1:2:w/2; elines=2:2:w/2;

    A = analysis_9_7(A',1,0)'; % horizontal analysis
    L = A(:,1:2:w); % horizontal L band

    % 3 lifting steps flipping high frequency spectra
    L(:,elines) = L(:,elines) + filtnflip(L(:,olines),M,Hb,Hl,beta);
    L(:,olines) = L(:,olines) - filtnflip(L(:,elines),M,Hb,Hl,beta);
    L(:,elines) = L(:,elines) + filtnflip(L(:,olines),M,Hb,Hl,beta);
    A(:,1:2:w) = L; % horizontal L band

    A = analysis_9_7(A,1,0); % vertical analysis

    % reverse lifting steps for the H band after DWT analysis
    HLolines = 1:4:w; HLelines = 3:4:w;
    A(2:2:h,HLelines) = A(2:2:h,HLelines) - A(2:2:h,HLolines);
    A(2:2:h,HLolines) = A(2:2:h,HLolines) + A(2:2:h,HLelines);
    A(2:2:h,HLelines) = A(2:2:h,HLelines) - A(2:2:h,HLolines);

    A = synthesis_9_7(synthesis_9_7(A,1,0)',1,0)';
end
if skipdwt==0, A = analysis_9_7(analysis_9_7(A',1,0)',1,0); end

Xout = Xin; Xout(hindices,windices)=A;