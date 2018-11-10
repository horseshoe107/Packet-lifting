function P = psnr(mse,range)
% P = psnr(mse,range*)
%
% Convert MSE values to PSNR. If not supplied, range defaults to 255.
if nargin==1
    range=255;
end
P = 10*log10(range^2 ./ mse);