function Y = zerodc(X,~)
% Y = zerodc(X,nyquist*)
%
% modify filters specified in X so that their dc component is zero.
% each column of X is treated as a separate filter to be corrected.
% if the dc gain of the filter is Gd
% when the filter length is odd, the centre tap is modified by -Gd
% when the filter length is even the 2 central filter taps are modified by
% -Gd/2
%
% if the nyquist flag is set, the filters will be modified so that they
% have both zero dc gain and zero nyquist gain.
% note that we can express the respective gains in terms of the sum of the
% cosets of samples, A = sum X(1:2:N) and B = sum X(2:2:N)
% Gd = A+B, Gn = A-B  ->  2A = Gd+Gn, 2B = Gd-Gn
% then to correct the filter, we modify each coset appropriately
N = size(X,1); % size of column vector
dcgain = sum(X,1);
if nargin==1
    if mod(N,2)==1 % if odd length, modify central tap
        X((N+1)/2,:) = X((N+1)/2,:) - dcgain;
    else % modify 2 most central taps
        X(N/2,:) = X(N/2,:) - dcgain/2;
        X(N/2+1,:) = X(N/2+1,:) - dcgain/2;
    end
else % set both dc and nyguist gain to zero
    acgain = sum(X(1:2:N),1) - sum(X(2:2:N),1);
    Aval  = (dcgain+acgain)/2; % subtract this from the sum of 1:2:N samples
    Bval  = (dcgain-acgain)/2; % subtract from the sum of 2:2:N samples
    switch mod(N,4)
        case 1
            X((N+1)/2,:) = X((N+1)/2,:) - Aval; % 1:2:N coset
            X((N+1)/2-1,:) = X((N+1)/2-1,:) - Bval/2; % 2:2:N coset
            X((N+1)/2+1,:) = X((N+1)/2+1,:) - Bval/2;
        case 3
            X((N+1)/2,:) = X((N+1)/2,:) - Bval; % 2:2:N coset
            X((N+1)/2-1,:) = X((N+1)/2-1,:) - Aval/2; % 1:2:N coset
            X((N+1)/2+1,:) = X((N+1)/2+1,:) - Aval/2;
        case 2
            X(N/2,:) = X(N/2,:) - Aval; % 1:2:N coset
            X(N/2+1,:) = X(N/2+1,:) - Bval;
        case 0
            X(N/2+1,:) = X(N/2+1,:) - Aval; % 1:2:N coset
            X(N/2,:) = X(N/2,:) - Bval;
    end
end
Y=X;