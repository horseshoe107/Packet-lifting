function [Lkernpar0,Hkernpar0] = inband_kernels(shifts,N,mode)
% [Lkernpar0,Hkernpar0] = inband_kernels(shifts,N*,mode*)
% 
% usage example: inband_kernels2([0:1/8:1/2])
%
% Computes (2N+1)-tap FIR shift kernels from approximations of cubic spline
% interpolation. The shifts required should be specified in the range
% [-0.5,0.5] - if larger shifts are required, this can be implemented by
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
b=2; % Bsh is a (2b+1)-tap filter
for m=1:length(shifts)
    sigma = shifts(m); % NB: require |sigma|<1
    for n=-b:b
        t = abs(n-sigma); % shift FORWARD means argument is (n-sigma)
        if t<1  Bsh(n+b+1) = 1/2*t^3 - t^2 + 2/3;
        else
            if t<2  Bsh(n+b+1) = -t^3/6 + t^2 - 2*t + 4/3;
            else    Bsh(n+b+1) = 0;
            end
        end
    end
    phiL(:,m) = conv(pL,Bsh);
    phiH(:,m) = conv(pH,Bsh);
end
C = (size(phiL,1)+1)/2; % central location of the L output

% construct kernels for (shift parity 0)
polyphiL = analysis_5_3(phiL,1,0);
polyphiH = analysis_5_3(phiH,1,0);
% subsample down
LtoL = polyphiL(C-Neven:2:C+Neven,:);
LtoH = polyphiL(C-Nodd:2:C+Nodd,:);
HtoL = polyphiH(C+1-Nodd:2:C+1+Nodd,:);
HtoH = polyphiH(C+1-Neven:2:C+1+Neven,:);
if normalisation% normalise polyphase components
    LtoL = LtoL ./ repmat(sum(LtoL,1),size(LtoL,1),1); % DC gain of 1
    LtoH(size(LtoH,1)/2:size(LtoH,1)/2+1,:) = ... % DC gain of 0, even length filter
        LtoH(size(LtoH,1)/2:size(LtoH,1)/2+1,:) - repmat(sum(LtoH,1)/2,2,1);
    HtoL(size(HtoL,1)/2:size(HtoL,1)/2+1,:) = ... % DC gain of 0, even length filter
        HtoL(size(HtoL,1)/2:size(HtoL,1)/2+1,:) - repmat(sum(HtoL,1)/2,2,1);
    nyqadj = ((-1).^(1:size(LtoL,1))*LtoL) ./ ((-1).^(1:size(HtoH,1))*HtoH);
    HtoH = HtoH .* repmat(nyqadj,size(HtoH,1),1); % same Nyquist gain as LtoL
end
% interleave to construct in-subband shift kernels
Lkernpar0(1+Npar:2:1+Npar+2*Neven,:) = LtoL;
Lkernpar0(2-Npar:2:2-Npar+2*Nodd,:)  = HtoL;
Hkernpar0(1+Npar:2:1+Npar+2*Neven,:) = HtoH;
Hkernpar0(2-Npar:2:2-Npar+2*Nodd,:)  = LtoH;
if any(mode=='W')
    fid = fopen('v2poly.dat','wt');
    fprintf(fid,'%d %d %d\n',size(Lkernpar0,1),numel(Lkernpar0)*4,2);
    writematrix(fid,Lkernpar0');
    writematrix(fid,Hkernpar0');
    fclose(fid);
    inband_kernels2(shifts,N,mode);
end