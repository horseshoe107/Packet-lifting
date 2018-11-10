function [G,S] = packlift_synthgains(synsteps,Ht,Hc,Hp)
% [G,S] = packlift_synthgains(synsteps,Ht*,Hc*,Hp*)
% 
% packlift_synthgains calculates the 1-dimensional synthesis vectors and
% the associated synthesis gains for excitations in packet lifting
% subbands. The sequence of dwt and packet lifting steps must be specified
% through the synsteps argument, which is a string composed of 'L', 'H',
% '2' and 'P' characters. Recognised control sub-strings are:
% 'L': 1 level of lowpass dwt 9x7 synthesis
% 'H': 1 level of highpass dwt 9x7 synthesis
% 'L2': synthesis from the LH' packet subband. Due to the nature of the
%       packet lifting structure, an excitation in LH' will stimulate
%       excitations in both LH and HL
% 'H2': synthesis from the HL' packet subband
% 'LP': synthesis from the LH' packet subband, with the addition that all 3
%       packet lifting steps are used (Hp is added)
% 'HP': synthesis from the HL' packet subband with all 3 packet lifting
%       steps
% these control sub-strings must be given in analysis order
if all(ismember(synsteps,'LH2P'))
    N = length(synsteps);
else error('Unknown control character entered')
end
if (nargin<2)
    Ht=1; end
if (nargin<3)
    Hc=fliplr(Ht); end
if (nargin<4)
    Hp=Hc; end
% for explanation of the offset requirement: (bk3,p7)
if mod(length(Ht),2)==0
    offset=1; % filters are odd order
else offset=-1;
end

n=N; S=1;
while n>0
    switch synsteps(n)
        case '2' % dwt; expand with zeros then convolve
            switch synsteps(n-1)
                case 'L' % if LH' with 2 lifting steps
                    [~,LHpath]=synthG('LH',Ht,Hc,0,S);
                    [~,HLpath]=synthG('HLT',Ht,Hc,0,S);
                case 'H' % if HL' with 2 lifting steps
                    [~,LHpath]=synthG('LHC',Ht,Hc,0,S);
                    [~,HLpath]=synthG('HL,HLTC',Ht,Hc,0,S);
            end
            S=sumvectors(LHpath,HLpath,offset);
            n=n-2;
        case 'P'
            switch synsteps(n-1)
                case 'L' % if LH' with precompensate step
                    [~,LHpath]=synthG('LH,LHPT',Ht,Hc,Hp,S);
                    [~,HLpath]=synthG('HLT',Ht,Hc,Hp,S);
                case 'H' % if HL' with precompensate step
                    [~,LHpath]=synthG('LHC,LHP,LHPTC',Ht,Hc,Hp,S);
                    [~,HLpath]=synthG('HL,HLTC',Ht,Hc,Hp,S);
            end
            S=sumvectors(LHpath,HLpath,offset);
            n=n-2;
        otherwise % expand through usual dwt vectors
            endn=n;
            while n>0 && not(ismember(synsteps(n),'2P'))
                n=n-1;
            end
            [~,S] = synthG(synsteps(n+1:endn),Ht,Hc,Hp,S);
    end
end
G = sum(S.^2);