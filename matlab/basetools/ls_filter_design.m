function B = ls_filter_design(N,hdhat,hdF,rhohat,rhoF,varargin)
% function B = ls_filter_design(N,hdhat,hdF,rhohat,rhoF,S*)
%
% Designs an FIR filter B of order N by minimising the weighted
% least-squares error. NOTE: N currently must be even.
% hdF and rhoF must be vectors of frequencies in the normalised range
% [0,1]. If hdhat consists of M values, then hdF must have length M+1.
% hdhat(n) specifies the desired, _constant_ response in the section
% [hdF(n),hdF(n+1)]. No polynomial or otherwise interpolation is used along
% the section. To define a smoother desired response, increase the length
% of both vectors.
% rhohat and rhoF must follow the same rules described for hd.
%
% Example of function use
% rhohat = [1 0 1]; % define weighting function
% rhoF = [0 1/3 2/3 1];
% hdhat = [1 0 0]; % define desired response
% hdF = [0 1/3 2/3 1];
% B = ls_filter_design(8,hdhat,hdF,rhohat,rhoF);

% setup
if mod(N,2)==1 % if N is odd
    error('Order N of filter must be even');
end
U = N/2; L = -U; % for odd order, U and L are non-integers
N = N+1; % from here on, N sets the number of taps
if nargin>=6
    S = varargin{1}; % total number of samples
else S = 6000; % default
end
deltaw = pi/S; % frequency spacing (deltaw->0 gives true integration)
w = deltaw * (0:(S-1)); % left hand indexing

% convert normalised frequencies to sample locations in w
hdF = round(hdF * S);
hdF = min(max(hdF,0),S); % shift inside range [0,S]
rhoF = round(rhoF * S);
rhoF = min(max(rhoF,0),S);

rhat = zeros(S,1); % by default, undefined regions have zero weight
for n = 1:length(rhohat)
    % each pass selects a new non-overlapping section
    rhat((rhoF(n)+1):rhoF(n+1)) = rhohat(n)^2;
end
dhat = rhat;
for n = 1:length(hdhat)
    dhat((hdF(n)+1):hdF(n+1)) = hdhat(n) * rhat((hdF(n)+1):hdF(n+1));
end

% numerical integration (check out inbuilt function ode45 later)
% exploit conjugate symmetry for integration over [0,pi]
d=zeros(N,1); r=zeros(N,1); R=zeros(N); % initialise
for n = L:U % inverse DTFT of d
    d(n-L+1) = real(exp(1j*w*n)*dhat) / S;
end % note e^jwn is a row vector, dhat is a column vector
for n = 0:(N-1) % inverse DTFT of r
    r(n+1) = real(exp(1j*w*n)*rhat) / S;
end
for n = 1:N % form R matrix from r
    for m = 1:N
        R(m,n) = r(abs(m-n)+1);
    end
end
B = (inv(R) * d)'; % LS-optimal coefficients, returned as a row vector