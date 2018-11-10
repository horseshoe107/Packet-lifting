x = 'P'; % '2' controls additional update/precomp step

% % Finding the 1st level packet lifting synth vectors is not particularly
% % easy. This is because the path is not a pure synthesis; the first step is
% % to analyse down to either LL/LH or HL/HH band pairs. The vectors produced
% % in these bands vary depending on the position of the unit excitation in
% % L or H (with periodicity 2)
% LtoL = [0.026749 -0.078223 0.602949 -0.078223 0.026749];
% LtoH = [0.045636 -0.295636 -0.295636 0.045636];
% HtoL = [-0.028772 0.557543 -0.028772];
% HtoH = [-0.016864 0.266864 0.266864 -0.016864];

% calculate 1-dimensional synthesis vectors
[G,L] = synthG('L');
[G,H] = synthG('H');
[G,LL] = synthG('LL');
[G,HL] = synthG('HL');
[G,LH] = synthG('LH');
[G,HH] = synthG('HH');
[G,LLL] = synthG('LLL');
[G,LLH] = synthG('LLH');
[G,HHL] = synthG('HHL');
[G,HHH] = synthG('HHH');

% 2nd level packet lifting synthesis vectors
[G,LH] = packlift_synthgains(['L' x],Ht,Hc,Hp);
[G,HL] = packlift_synthgains(['H' x],Ht,Hc,Hp);
% 3rd level packet lifting synthesis vectors
[G,LHL] = packlift_synthgains(['L' x 'L'],Ht,Hc,Hp);
[G,LHH] = packlift_synthgains(['L' x 'H'],Ht,Hc,Hp);
[G,HLL] = packlift_synthgains(['H' x 'L'],Ht,Hc,Hp);
[G,HLH] = packlift_synthgains(['H' x 'H'],Ht,Hc,Hp);

% Find largest correlation between synthesis vectors for varying shifts
clear Y, Y{1}=LL; Y{2}=HL; Y{3}=LH; Y{4}=HH;
A = zeros(4);
% Both LH' and HL' bands end up centred around the same dominant
if x == 'P' % path (bk3,p10-11)
    dompath = 3; % LH dominates
else % ''/'U': ordinary 2 step lifting, or additional Update
    dompath = 2; % HL dominates
end
for n=1:4
    for m=1:4
        if (n==2||n==3) a=dompath; else a=n; end
        if (m==2||m==3) b=dompath; else b=m; end
        T = translates(Y{n},Y{m},4,b-a);
        if max(T)<max(abs(T))
            A(n,m) = min(T);
        else A(n,m) = max(T);
        end
    end
end
% Test for deep packet transform
Y{1}=LLL; Y{2}=HLL; Y{3}=LHL; Y{4}=HHL;
Y{5}=LLH; Y{6}=HLH; Y{7}=LHH; Y{8}=HHH;
B = zeros(8);
for n=1:8
    for m=1:8
        if (n==2)|(n==3) a=dompath; else a=n; end
        if (n==6)|(n==7) a=dompath+4; else a=n; end
        if (m==2)|(m==3) b=dompath; else b=m; end
        if (m==6)|(m==7) b=dompath+4; else b=m; end
        T = translates(Y{n},Y{m},8,b-a);
        if max(T)<max(abs(T))
            B(n,m) = min(T);
        else B(n,m) = max(T);
        end
    end
end

% calculate correlation coefficients
C = zeros(4);
for n=1:4
    for m=1:4
        C(n,m) = A(n,m)/sqrt(A(n,n)*A(m,m));
    end
end
D = zeros(8);
for n=1:8
    for m=1:8
        D(n,m) = B(n,m)/sqrt(B(n,n)*B(m,m));
    end
end
disp(C), disp(D)