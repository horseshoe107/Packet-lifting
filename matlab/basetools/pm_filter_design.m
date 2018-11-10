function [varargout] = pm_filter_design(F,W,stopripple,varargin)
% [H] = pm_filter_design(F,W,stopripple)
% [Ffor,Frev] = pm_filter_design(F,W,stopripple,mode*,thresh*)
%                             * -  optional arguments
%
% Designs a linear phase FIR LPF `H' with even order (odd number of taps).
% This constraint was introduced in order to simplify filtering.
% The order is progressively increased until the filter meets the stopband
% ripple requirement (and the associated passband ripple requirement, set
% by the weighting argument W)
% pm_filter_design will additionally factorise `H' by min/max phase
% components if 2 outputs are required.
%
% Note that the passband is designed to have unity gain, which means that
% in general the filter will not have DC gain 1.
%
% The mode argument is an optional string. If the char 'F' is asserted, the
% filter will be designed to be factorisable: that is, zeros that would
% normally appear in conjugate pairs at the stopband will appear in
% quadruples. This is done by forcing the designed response to be strictly
% positive.
% If '-' is asserted, the response will be forced to be negative in the
% stopband region. There will be one lone zero pair on the unit circle.
%
% Zeros assumed to correspond to the stopband are moved to the unit circle
% after factorisation. By default, this occurs if the original magnitude of
% the zero is within the range [0.95,1.05]. The 5th argument will instead
% set this threshold to [1-thresh,1+thresh]

A = zeros(1,length(F));
if nargin>=4
    mode = varargin{1};
    if any(mode=='F'), A=A+stopripple; end
    if any(mode=='-'), A=A-stopripple; end
else mode='';
end
if nargin==5
    thresh = varargin{2};
else thresh = 0.05;
end
A(1:2)=[1 1];

N=4;
while 1
    N = N+2;
    [H,~,res] = firpm(N,F,A,W);
    errveclen = length(res.error);
    % Assume the 3rd freq value is start of the stopband specification
    stoperrors = (floor(F(3)*errveclen):errveclen);
    if max(abs(res.error(stoperrors))) <= stopripple % (stopripple*sum(H))
        break
    end
end
ERR = res.error;

if nargout==1
    varargout{1} = H;
end
if nargout==2
    % set boundaries for what is considered a passband/stopband zero
    up = 1+thresh; low=1-thresh;
    if (any(mode=='F') || any(mode=='-'))
        E = sqrt(1+max(ERR));
        R = roots(H); Rmag=abs(R);
        I1pass = find(Rmag>=up); I1stop = find((up>Rmag)&(Rmag>=1));
        I2pass = find(low>=Rmag); I2stop = find((1>Rmag)&(Rmag>low));
    else
        E = sqrt(1+max(ERR));
        R = roots(H); Rmag=abs(R);
        I1pass = find(Rmag>=up); I2pass = find(low>=Rmag);
        Istop = find((up>Rmag)&(Rmag>low));
        I1stop = Istop([1:4:length(Istop) 2:4:length(Istop)]);
        I2stop = Istop([3:4:length(Istop) 4:4:length(Istop)]);
    end
    A = poly([R(I1pass) ; exp(1i*angle(R(I1stop)))]); A=real(A); Ffor=A/sum(A)*E;
    B = poly([R(I2pass) ; exp(1i*angle(R(I2stop)))]); B=real(B); Frev=B/sum(B)*E;
    varargout{1} = Ffor; varargout{2} = Frev;
end