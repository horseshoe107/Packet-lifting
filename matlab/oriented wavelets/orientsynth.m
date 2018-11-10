function orientsynth(mode)
% orientsynth(mode)
%
% Calculate synthesis gains and correlation coefficients between synthesis
% vectors produced by oriented variants of the 5/3 dwt.
%
% Mode 0: Gains from a 1D oriented wavelet synthesis
% Mode 1: Gains from 2D - ordinary dwt, followed by oriented synthesis
% Mode 2: Correlation coefficients between 2D synthesis patterns
switch mode
    case 0
        disp('1D oriented synthesis gains')
        disp('L           H')
    case 1
        disp('2D oriented synthesis gains')
        disp('1st letter: normal horizontal dwt,  2nd letter: oriented vertical dwt')
        disp('    LL        LH        HL        HH')
    case 2
        disp('Correlation coefficients between pairs of 2D synthesis functions')
        disp('1st letter: normal horizontal dwt,  2nd letter: oriented vertical dwt')
        disp('    LL-LH     LL-HL     LL-HH     LH-HL     LH-HH     HL-HH')
end

% DATA
% L synthesis vector (row vector)
Lhor = [0.5000 1.0000 0.5000];
% H synthesis vector
Hhor = [-0.2500 -0.5000 1.5000 -0.5000 -0.2500];
% interpolation kernels for varying degrees of shift
K_LUT = [ 0.0000000000 -0.0000000000  0.0000000000  1.0000000000  0.0000000000 -0.0000000000  0.0000000000 
-0.0057314933  0.0213902240 -0.0798294029  0.9674625206  0.1198302579 -0.0315853847  0.0084632783 
-0.0088233996  0.0329293754 -0.1228941022  0.8796881878  0.2687588056 -0.0678352739  0.0181764069 
-0.0097770136  0.0364883116 -0.1361762329  0.7514540266  0.4325511590 -0.1018238761  0.0272836254 
-0.0090991941  0.0339586548 -0.1267354249  0.5973263672  0.5973263672 -0.1267354249  0.0339586548 ];

for shift=0:1/8:1/2
    K = K_LUT(1+shift*8,:);
    % effective filters for repeated shift (k*k) and return (k*~k)
    R = conv(K,fliplr(K)); % return filter
    D = conv(K,K);         % double shift filter

    % refer to bk3p53
    % L synthesis rows
    L0 = 1; % centre row
    L1 = K/2; % -1, 1 rows
    % H synthesis rows
    H0 = -0.5*R;
    centretmp = (size(R,2)+1)/2;
    H0(:,centretmp) = H0(:,centretmp)+2; % centre row
    H1 = -K/2;
    H2 = -D/4;
    if mode==0
        Lsynthgain1D = sum(L0.^2,2) + 2*sum(L1.^2,2);
        Hsynthgain1D = sum(H0.^2,2) + 2*sum(H1.^2,2) + 2*sum(H2.^2,2);
        disp([Lsynthgain1D Hsynthgain1D])
    end
    % 2D synthesis impulse response rows
    LL0 = conv2(Lhor,L0); % 1st letter denotes horizontal, normal dwt synthesis
    LL1 = conv2(Lhor,L1); % 2nd letter denotes vertical, oriented dwt
    HL0 = conv2(Hhor,L0); % numeral denotes row (0 is centre row)
    HL1 = conv2(Hhor,L1);
    LH0 = conv2(Lhor,H0);
    LH1 = conv2(Lhor,H1);
    LH2 = conv2(Lhor,H2);
    HH0 = conv2(Hhor,H0);
    HH1 = conv2(Hhor,H1);
    HH2 = conv2(Hhor,H2);
    LLgain = sum(LL0.^2,2) + 2*sum(LL1.^2,2);
    HLgain = sum(HL0.^2,2) + 2*sum(HL1.^2,2);
    LHgain = sum(LH0.^2,2) + 2*sum(LH1.^2,2) + 2*sum(LH2.^2,2);
    HHgain = sum(HH0.^2,2) + 2*sum(HH1.^2,2) + 2*sum(HH2.^2,2);
    if mode==1
        disp([LLgain LHgain HLgain HHgain])
    end
    if mode==2
        % correlation between L-Lorient and L-Horient
        dotLL_LH = translates(LL1,LH2) + ...
            translates(LL0,LH1) + ...
            translates(fliplr(LL1),LH0);
        dotLL_HL = translates(LL1,HL1,0,1) + ...
            translates(LL0,HL0,0,1) + ...
            translates(fliplr(LL1),fliplr(HL1),0,1);
        dotLL_HH = translates(LL1,HH2,0,1) + ...
            translates(LL0,HH1,0,1) + ...
            translates(fliplr(LL1),HH0,0,1);
        dotLH_HL = translates(LH2,HL1,0,1) + ...
            translates(LH1,HL0,0,1) + ...
            translates(LH0,fliplr(HL1),0,1);
        dotLH_HH = translates(LH2,HH2,0,1) + translates(LH1,HH1,0,1) + ...
            translates(LH0,HH0,0,1) + translates(fliplr(LH1),fliplr(HH1),0,1) + ...
            translates(fliplr(LH2),fliplr(HH2),0,1);
        dotHL_HH = translates(HL1,HH2) + ...
            translates(HL0,HH1) + ...
            translates(fliplr(HL1),HH0);
        corrLL_LH = dotLL_LH / sqrt(LLgain*LHgain);
        corrLL_HL = dotLL_HL / sqrt(LLgain*HLgain);
        corrLL_HH = dotLL_HH / sqrt(LLgain*HHgain);
        corrLH_HL = dotLH_HL / sqrt(LHgain*HLgain);
        corrLH_HH = dotLH_HH / sqrt(LHgain*HHgain);
        corrHL_HH = dotHL_HH / sqrt(HLgain*HHgain);
        disp([corrLL_LH corrLL_HL corrLL_HH corrLH_HL corrLH_HH corrHL_HH])
    end
end