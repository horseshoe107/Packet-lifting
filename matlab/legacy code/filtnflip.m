function Y = filtnflip(X,M,Hb,Hl,beta)
% function Y = filtnflip(X,M,Hb,Hl,beta)
%
% function to select out a desired portion of spectrum, and flip that about
% its own centre. (ie, not generally pi/2)
% the steps in order are:
% 1. Upsample by M
% 2. Very fine filtering of a flipped image by band pass filter Hb
% 3. Modulation by cos(beta*pi*n)
% 4. Lazy LPF Hl keeping only the baseband component
% 5. Downsample by M
%
% The size of M determines a tradeoff between the complexity of Hb and Hl.
% M can be as low as 3. In general, M should be an odd number. Picking an
% even number increases the complexity of Hb without a corresponding saving
% in Hl.

% HbN=Hb(1,:); HlN=Hl(1,:); % extract filter coefficients
% if (size(Hb,1)>=2) HbD=Hb(2,:); else HbD=1; end % FIR/IIR
% if (size(Hl,1)>=2) HlD=Hl(2,:); else HlD=1; end

[h,w] = size(X);
% upsample by M
U = zeros(h*M,w); U(1:M:h*M,:)=X;
% bandpass filter
Xb = filtextend(U,Hb,1,1,[1,h*M]);
% modulate to baseband
modker = cos(pi*beta*(1:h*M))';
Xm = Xb .* repmat(modker,1,w);
% low pass filter
Xl = filtextend(Xm,Hl,1,2,[1,h*M])*2; % restore after filtering half of the modulated signal
% downsample
D = Xl(1:M:h*M,:)*M; % restore for upsampling factor
Y = D;

% % diagnostics
% % shows effect of band pass filter - are the dc spikes rejected?
% figure, semilogy([abs(fftpsd(U)) abs(fftpsd(Xb))])
% % effect of low pass filter - reject the upshifted copy
% figure, semilogy([abs(fftpsd(Xm)) abs(fftpsd(Xl))])