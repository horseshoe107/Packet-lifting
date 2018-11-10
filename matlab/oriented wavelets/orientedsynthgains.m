function orientedsynthgains(shifts)
% orientedsynthgains(shifts)
%
% Generate oriented synthesis vectors from an arbitrary orientation field,
% described by a set of relative shifts between adjacent rows. The 5/3 dwt
% is used as the base for generation.
% The shifts argument should be a 4-element vector, where the order of the
% elements corresponds to adjacent row-pairs as they increase. The L
% excitation is centred on row 2, while H excitation is centred on row 3
%
% ------ 1st row
%                <- 1st shift
% ------ 2nd row
%                <- 2nd shift
% ------ 3rd row     etc
% ------ 4th row
% ------ 5th row

Lhor = [0.5000 1.0000 0.5000]; % L synthesis vector (row vector)
Hhor = [-0.2500 -0.5000 1.5000 -0.5000 -0.2500]; % H synthesis vector

% refer to bk3p53
% L synthesis rows
L0 = 1; % centre row
L1 = fetchkernel(shifts(1))/2;
L_1 = fetchkernel(-shifts(2))/2;
% H synthesis rows
H1 = -fetchkernel(shifts(2))/2;
Dtmp = conv( fetchkernel(shifts(2)) , fetchkernel(shifts(1)) );
H2 = -Dtmp/4;
R1 = conv(fetchkernel(shifts(2)),fetchkernel(-shifts(2)));
R2 = conv(fetchkernel(shifts(3)),fetchkernel(-shifts(3)));
H0 = -(R1+R2)/4;
H0((size(H0,2)+1)/2) = H0((size(H0,2)+1)/2)+2; % centre sample
H_1 = -fetchkernel(-shifts(3))/2;
Dtmp = conv( fetchkernel(-shifts(3)) , fetchkernel(-shifts(4)) );
H_2 = -Dtmp/4;

L = stackrows(L1,L0,L_1);
H = stackrows(H2,H1,H0,H_1,H_2);
disp(L)
disp(H)
Lsynthgain1D = sum(sum(L.^2));
Hsynthgain1D = sum(sum(H.^2));
disp('1D synthesis gains: L, H');
disp([Lsynthgain1D Hsynthgain1D])

% 2D synth gain calculations
LL = conv2(Lhor,L);
HL = conv2(Hhor,L);
LH = conv2(Lhor,H);
HH = conv2(Hhor,H);
LLsynthgain = sum(sum(LL.^2));
HLsynthgain = sum(sum(HL.^2));
LHsynthgain = sum(sum(LH.^2));
HHsynthgain = sum(sum(HH.^2));
disp('2D synthesis gains:');
disp('L dwt, L oriented - L dwt, H oriented - H dwt, L oriented - H dwt, H oriented');
disp([LLsynthgain LHsynthgain HLsynthgain HHsynthgain])