function [LLkern,H1kern,LHkern,H2kern] = packetband_kernels(shifts,N,mode)
% [LLkern,H1kern,LHkern,H2kern] = packetband_kernels(shifts,N*,mode*)
% 
% Computes (2N+1)-tap FIR shift kernels from approximations of cubic spline
% interpolation. The shifts required should be specified in the range
% [z,z+0.5] where z is any integer; any shifts in [z-0.5,z] can be
% implemented by the mirror image filter of the corresponding positive
% shift.
%
% The shift kernels computed act directly on subband coefficients that are
% 2 dwt levels deep; ie, they operate on the LL, H and LH bands. For any
% specific fractional-pixel shift, there are 16 possible kernels. There are
% 4 possible subbands when writing out to the destination pixel. There are
% also 4 different subband contexts the filter will act on, depending on
% the integer value of the shift.
%
% The mode argument is a string, with chars flagging desired modes. 'N'
% will activate normalisation for all the polyphase filters. This mode is
% active by default, so a '' argument will disable this.
%
% See bk3 p22-23 for a description of the shift kernel construction
if nargin<2, N=10; end % length of the FIR shift kernels are set to 2N+1 taps
if nargin<3, mode='N'; end
normalisation=any(mode=='N');

z = round(mean(shifts)); % *guess* whole-pixel shift
if any(abs(shifts-z)> 0.5)
    error('shifts must be in the range [z-0.5,z+0.5] where z is any integer');
end
if z<0
    error('shifts must be greater than -0.5 (integer component must be non-negative)');
end
shifts = shifts - z;

% set up impulse excitations, then synthesise to image domain
C = 33; % central location of the LL excitation
for i=1:4
    imp{i}=zeros(2*C-1,1); imp{i}(C+i-1,1)=1;
    % LL->33, H1-34, LH-35, H2-36
    % Note, a H excitation will not be affected by the inner synthesis
    synth{i}=synthesis_5_3(synthesis_5_3(imp{i},2,0),1,0);
end
C = C+z; % overwritten: central location of the LL *output*

% apply cubic b-spline shift: first, calculate the b-spline coefficients
% by causal and anticausal IIR filtering
p1 = -2 + sqrt(3); % causal pole
a0 = (1-p1); % scaling factor so filter has dc gain 1
for i=1:4
    p{i} = filter(a0,[1 -p1],synth{i});
    p{i} = flipud(filter(a0,[1 -p1],flipud(p{i})));
end

% calculate bspline values at shifted pixel positions
b=2; % Bspline has a nonzero support of 2b+1 taps
for m=1:length(shifts)
    sigma = shifts(m); % NB: require |sigma|<1
    for n=(-b:b)
        t = abs(n-sigma); % shift FORWARD means argument is (n-sigma)
        if t<1  Bsh(n+b+z+1) = 1/2*t^3 - t^2 + 2/3;
        else
            if t<2  Bsh(n+b+z+1) = -t^3/6 + t^2 - 2*t + 4/3;
            else    Bsh(n+b+z+1) = 0;
            end
        end
    end
    for i=1:4
        tmp = conv(p{i},Bsh);
        phi{i}(:,m) = tmp((1+b):length(tmp)); % drop the first b samples
        polyphi{i}  = analysis_5_3(analysis_5_3(phi{i},1,0),2,0);
    end
end


if normalisation
    for m=1:length(shifts)
        % set DC gain to 1 for LL->LL cofilter
        polyphi{1}(C-N+mod(N-z,4):4:C+N,m) = polyphi{1}(C-N+mod(N-z,4):4:C+N,m)...
            / sum(polyphi{1}(C-N+mod(N-z,4):4:C+N,m));
        % set DC and Nyg gain to 0 for LL->H cofilter
        tmpindex = C-N+mod(N+1-z,2):2:C+N;
        polyphi{1}(tmpindex,m) = zerodc(polyphi{1}(tmpindex,m),'dc+ac=0');
        % set DC gain 0 for LL->LH cofilter
        tmpindex = C-N+mod(N+2-z,4):4:C+N;
        polyphi{1}(tmpindex,m) = zerodc(polyphi{1}(tmpindex,m));
        % set DC gain 0 for H->LL cofilter
        tmpindex = C+1-N+mod(N-z-1,4):4:C+1+N;
        polyphi{2}(tmpindex,m) = zerodc(polyphi{2}(tmpindex,m));
        tmpindex = C+3-N+mod(N-z-3,4):4:C+3+N;
        polyphi{4}(tmpindex,m) = zerodc(polyphi{4}(tmpindex,m));
        % set DC gain 0 for LH->LL cofilter
        tmpindex = C+2-N+mod(N-z-2,4):4:C+2+N;
        polyphi{3}(tmpindex,m) = zerodc(polyphi{3}(tmpindex,m));
        % set DC gain 0 for LH->H cofilter
        tmpindex = C+2-N+mod(N-z+1-2,2):2:C+2+N;
        polyphi{3}(tmpindex,m) = zerodc(polyphi{3}(tmpindex,m));
    end
end

% interleave to construct in-subband shift kernels
for i=0:3 % kernels organised according to *destination* pixel
    j=mod(-i-z,4); % (0,3,2,1 for z=0), (3,2,1,0 for z=1)
    k=mod(-i-z+1,4); % 1,0,3,2
    l=mod(-i-z+2,4); m=mod(-i-z+3,4);
    LLkern((mod(N+i,4)+1):4:(2*N+1),:) = polyphi{j+1}(C+j-N+mod(N+i,4):4:C+j+N,:);
    H1kern((mod(N+i,4)+1):4:(2*N+1),:) = polyphi{k+1}(C+k-N+mod(N+i,4):4:C+k+N,:);
    LHkern((mod(N+i,4)+1):4:(2*N+1),:) = polyphi{l+1}(C+l-N+mod(N+i,4):4:C+l+N,:);
    H2kern((mod(N+i,4)+1):4:(2*N+1),:) = polyphi{m+1}(C+m-N+mod(N+i,4):4:C+m+N,:);
end

if any(mode=='W')
    if z==0 % create fresh file, and add header
        fid = fopen('v4pack.dat','wt');
        fprintf(fid,'%d %d %d\n',size(LLkern,1),numel(LLkern)*16,4);
    else % append
        fid = fopen('v4pack.dat','at');
    end
    writematrix(fid,LLkern');
    writematrix(fid,H1kern');
    writematrix(fid,LHkern');
    writematrix(fid,H2kern');
    fclose(fid);
    if z~=3
        packetband_kernels(shifts+z+1,N,mode);
    end
end