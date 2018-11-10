function [out,T1,T2,T3] = symbolic(f,sigma,N)
% [out,T1,T2,T3] = symbolic(f,sigma,N*)
%
% This function takes a single normalised (in range [0,1]) frequency f and
% a shift sigma - which denotes the amount of shift existing between
% adjacent rows on the hypothetical test image.
% The frequency component f in the low-pass subband after the transform is
% performed is composed of a scaled version of x(pi*w) and a scaled version
% of x(pi*w - pi) where x is the input image. These 2 scaling coefficients
% are returned in the row vector out.
% T1-T3 are the intermediate values of the scaling coefficients after each
% lifting step in the transform.
%
% Refer to p84-85 of book 3 for further notes

breakpoints = [0.375 0.4375 0.5625 0.625];
% breakpoints = [0.425 0.5 0.625 0.7];
if (f<0)||(f>1)
    error('Frequency must be in the normalised range [0,1]');
end
w = f*pi;

h0_53 = [-1/8 1/4 3/4 1/4 -1/8];   h0_N = -2:2;
g0_53 = [1/2 1 1/2];               g0_N = -1:1;

S = diag([exp(-1j*sigma*w) exp(-1j*sigma*(w-pi))]); % 2x2 shift matrix
if nargin==3
    B = diag([bpflookup(breakpoints,f,N) bpflookup(breakpoints,f-1,N)]);
else
    B = diag([bpflookup(breakpoints,f) bpflookup(breakpoints,f-1)]);
end
H = 0.5*[h0_53*exp(-1j*h0_N'*w) h0_53*exp(-1j*h0_N'*(w-pi))]; % 1x2 analysis matrix
G = [g0_53*exp(-1j*g0_N'*w) g0_53*exp(-1j*g0_N'*(w-pi))]'; %2x1 synth matrix
I = eye(2); % 2x2 identity matrix

% 1st lifting step - lowpass rows have horizontal spectra removed
% x'_2k = T1 * x_2k
% x'_2k1 = x_2k1 (odd rows unchanged)
T1 = I - G*H*B;

% 2nd lifting step
% x''_2k = x'_2k (even rows unchanged)
% x''_2k1 = T2 * x_2k
% x''_2k_1 = T2_ * x_2k
T2 = S*(I - 0.5*T1 - 0.5*(S^-2)*T1*S^2);
T2_ = (S^-1)*(I - 0.5*T1 - 0.5*S^2*T1*(S^-2));

% 3rd lifting step
% x'''_2k1 = x''_2k1 (odd rows unchanged)
% x'''_2k = T3 * x_2k
T3 = T1*(I + 0.25*S*T2_ + 0.25*(S^-1)*T2);

% T3 is a 2x2 matrix that shows how our structure maps from a pair of
% points (the signal x at frequencies w and w-pi) to another pair of points
% (the baseband output, again at frequencies w and w-pi)
% To observe the subband output, we must multiply by the analysis matrix H,
% yielding a 1x2 transform
out = H*T3; % first component is signal, 2nd is aliasing