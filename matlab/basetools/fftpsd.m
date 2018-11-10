function [Y] = fftpsd(X)
% Use the fft to find the approximate power spectral density of X. X is
% windowed with the raised cosine before applying the fft. The returned
% vector Y is the squared magnitudes of the fft coefficients in [0,pi].
% (Note that Y is a symmetrical funcion if X is real valued)
%
% If X is a matrix, the fft is applied one-dimensionally to the columns and
% we average across the columns for Y.
% columns.
[h w] = size(X);
if (h==1)&&(w~=1)
  X = X';
  h = w;
  w = 1;
end
hsym = (h-1)/2;
coswin1D = cos(pi*(-hsym:hsym)/hsym)';
coswin = repmat((coswin1D+1)/2,1,w);
T = fft(coswin .* X);
Y = sum(abs(T).^2,2)/(w/2);
Y = Y(1:ceil(h/2),:);
if nargout==0
    figure, semilogy(Y);
    clear Y;
end