function [Lkernpar1,Hkernpar1] = inband_kernels2(shifts,N,mode)
% [Lkernpar1,Hkernpar1] = inband_kernels2(shifts,N*,mode*)
%
% usage example: inband_kernels2([1:1/8:3/2]) <- there may be an error here
% 
% Computes (2N+1)-tap FIR shift kernels from approximations of cubic spline
% interpolation. The shifts required should be specified in the range
% [1,1.5] - if larger shifts are required, this can be implemented by
% moving the sliding window filter to a different location to absorb any
% integer shift component.
%
% The shift kernels computed act directly on subband coefficients -
% however, because the subband coefficients are polyphase with periodicity
% 2, there are 4 possible filters applicable for each shift. At
% implementation, filter selection is based on:
% 1. Whether the destination pixel is in the L or H subband (output of
% interest is Lkern*** or Hkern***)
% 2. Whether the shift to be implemented, when rounded to the nearest
% integer, is odd or even. Within this function this is called the shift
% parity (output of interest is ***par0 or ***par1)
%
% The mode argument is a string, with chars flagging desired modes. 'N'
% will activate normalisation for all the polyphase filters. This mode is
% active by default, so a '' argument will disable this.
%
% See bk3 p22-23 for a description of the shift kernel construction
if nargin<2, N=4; end % length of the FIR shift kernels are set to 2N+1 taps
if nargin<3, mode='N'; end
normalisation=any(mode=='N');
if mod(N,2)==0
    Neven=N; Nodd=N-1; Npar=0;
else Neven=N-1; Nodd=N; Npar=1;
end

% set up impulse excitations in the L and H subbands
impL=zeros(61,1); impL(31,1)=1; % NB: vectors *must* be odd length
impH=zeros(61,1); impH(32,1)=1;
% synthesise to image domain
imL = synthesis_5_3(impL,1,0);
imH = synthesis_5_3(impH,1,0);

% apply cubic b-spline shift: first, calculate the b-spline coefficients
% by causal and anticausal IIR filtering
p1 = -2 + sqrt(3); % causal pole
a0 = (1-p1); % scaling factor so filter has dc gain 1
pL = filter(a0,[1 -p1],imL);
pL = flipud(filter(a0,[1 -p1],flipud(pL)));
pH = filter(a0,[1 -p1],imH);
pH = flipud(filter(a0,[1 -p1],flipud(pH)));

% calculate bspline values at shifted pixel positions
b=2; % Bspline has a nonzero support of 2b+1 taps
for m=1:length(shifts)
    sigma = shifts(m); % NB: require |sigma|<1
    for n=(-b:b)
        t = abs(n-sigma); % shift FORWARD means argument is (n-sigma)
        if t<1  Bsh(n+b+2) = 1/2*t^3 - t^2 + 2/3;
        else
            if t<2  Bsh(n+b+2) = -t^3/6 + t^2 - 2*t + 4/3;
            else    Bsh(n+b+2) = 0;
            end
        end
    end
    phiL(:,m) = conv(pL,Bsh);
    phiH(:,m) = conv(pH,Bsh);
end
C = (size(phiL,1)/2)+1; % central location of the L output
polyphiL = analysis_5_3(phiL,1,0);
polyphiH = analysis_5_3(phiH,1,0);

% construct kernels for (shift parity 1)
LtoL = polyphiL(C-Nodd:2:C+Nodd,:);
LtoH = polyphiL(C-Neven:2:C+Neven,:);
HtoL = polyphiH(C+1-Neven:2:C+1+Neven,:);
HtoH = polyphiH(C+1-Nodd:2:C+1+Nodd,:);
if normalisation % normalise
    LtoL = LtoL ./ repmat(sum(LtoL,1),size(LtoL,1),1); % DC gain of 1
    LtoH(ceil(size(LtoH,1)/2),:) = ... % DC gain of 0, odd length filter
        LtoH(ceil(size(LtoH,1)/2),:) - sum(LtoH,1);
    HtoL(ceil(size(HtoL,1)/2),:) = ... % DC gain of 0, odd length filter
        HtoL(ceil(size(HtoL,1)/2),:) - sum(HtoL,1);
    % Doesn't work because of divide by zero; stuff it
    % nyqadj = ((-1).^(1:size(LtoL,1))*LtoL) ./ ((-1).^(1:size(HtoH,1))*HtoH);
    % HtoH = HtoH .* repmat(nyqadj,size(HtoH,1),1); % same Nyquist gain as LtoL
end
% interleave to construct in-subband shift kernels
Lkernpar1(1+Npar:2:1+Npar+2*Neven,:) = HtoL;
Lkernpar1(2-Npar:2:2-Npar+2*Nodd,:)  = LtoL;
Hkernpar1(1+Npar:2:1+Npar+2*Neven,:) = LtoH;
Hkernpar1(2-Npar:2:2-Npar+2*Nodd,:)  = HtoH;
if any(mode=='W')
    fid = fopen('v2poly.dat','at');
    writematrix(fid,Lkernpar1');
    writematrix(fid,Hkernpar1');
    fclose(fid);
end