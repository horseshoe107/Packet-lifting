function skern = bspline_kernels(shifts,N,mode)
% skern = bsplines(shifts,N*,mode*)
% 
% Computes (2N+1)-tap FIR shift kernels from approximations of cubic spline
% interpolation. The shifts required should be in the range [-0.5,0.5] as
% the kernels become increasingly inaccurate for larger shifts.
%
% All splines are even symmetrical, ie Bk(t) = Bk(-t)
% B0(t) = 1                            in [0,1/2]
% B1(t) = -t+1                         in [0,1]
% B2(t) = -t^2 + 3/4                   in [0,1/2]
%       = 0.5*t^2 - 3/2*t + 9/8        in [1/2,3/2]
% B3(t) = 1/2*t^3 - t^2 + 2/3          in [0,1]
%       = -t^3/6 + t^2 - 2*t + 4/3     in [1,2]
%
% Note that for just discrete-time samples, B3[n] = 1/6,2/3,1/6
% To find the roots, we are effectively solving the quadratic 1+4x+x^2=0,
% which is -b +- sqrt(b^2 - 4ac)
%          --------------------- = -2 +- sqrt(3)
%                  2a
%
% All kernels are automatically normalised for DC gain 1

% default length of the FIR shift kernels if N not supplied
if nargin<2, N=3; end
if nargin<3, mode=''; end

p1 = -2 + sqrt(3); % causal pole
a0 = (1-p1); % scaling factor so filter has dc gain 1
imp=zeros(1,25); imp(1)=1; % set up arbitrary length impulse function
p = filter(a0,[1 -p1],imp); % filter for impulse response
p = conv(p,fliplr(p))/sum(p); % convolve with anticausal impulse response

b=2;
% calculate bspline values at shifted pixel positions
for m=1:length(shifts)
    sigma = shifts(m);
    for n=-b:b
        t = abs(n-sigma); % shift FORWARD means argument is (n-sigma)
        if t<1  Bsh(n+b+1) = 1/2*t^3 - t^2 + 2/3;
        else
            if t<2  Bsh(n+b+1) = -t^3/6 + t^2 - 2*t + 4/3;
            else    Bsh(n+b+1) = 0;
            end
        end
    end
    phi = conv(p,Bsh); % equivalent interpolation kernel
    C = (length(phi)+1)/2; % location of centre
    skern(:,m) = phi(C-N:C+N) / sum(phi(C-N:C+N)); % N-tap normalisation
end

if any(mode=='W')
    fid = fopen('baseker.dat','wt');
    fprintf(fid,'%d %d %d\n',size(skern,1),numel(skern),1);
    writematrix(fid,skern');
    fclose(fid);
end