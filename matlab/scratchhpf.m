% hpf_prelift adaptivity cases
% Case 1: 
s = 0.01;
w = 0:s:pi;
delta_shift = 1/16;
h0 = 3/4 + cos(w)/2 - cos(2*w)/4; % g0 = 1 + cos(w);
hpf = 38/64 - 19/32*cos(w) -5/32*cos(2*w) + 3/32*cos(3*w) +...
    4/32*cos(4*w) -2/32*cos(6*w);
x = sqrt(atan(pi./w)./w); x(1)=x(2); % can't let x(1) be Inf
S = 2*h0.*hpf.*x;
D = h0.*hpf.*x.*w*delta_shift;
sum_mag = sum(S)*s;
diff_mag = sum(D)*s;

% hpf_update adaptivity cases
s = 0.01;
w = 0:s:pi;
x = atan(pi./w)./w; x(1)=x(2); % can't let x(1) be Inf
x_predict_error = atan(pi./w).*w; % *2*delta_shift, but constant scale factor doesn't matter
h0 = 3/4 + cos(w)/2 - cos(2*w)/4; % g0 = 1 + cos(w);
hpf = 38/64 - 19/32*cos(w) -5/32*cos(2*w) + 3/32*cos(3*w) +...
    4/32*cos(4*w) -2/32*cos(6*w);
A = x_predict_error.*h0.^2;
B = x_predict_error.*h0.^2.*hpf.^2;
L_mean = sum(A)*s;
Alias_mean = sum(B)*s; % ratio of roughly 2x, with random data
C = x.*h0.^2.*hpf.^2;
D = x.*h0.^2.*hpf.^4;
figure(1), plot(w,[A' B' C' D'])


rate = [0.1 0.2 0.4 0.6 0.8 1];
% barbara
% Ordinary dwt MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A1 = psnr([238.976 141.517 67.0499 37.6582 24.508 16.6328])';
% Packet lifting MSEs adaptive w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A2 = psnr([244.532 147.277 75.059 43.8297 29.6309 20.9183])';
% Oriented wavelet MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A3 = psnr([227.144 125.421 55.764 31.1543 20.8333 13.6615])';
% 2 level oriented MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A4 = psnr([226.534 125.249 56.6402 31.455 20.8683 13.8522])';
% Oriented packet wavelet MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A5 = psnr([208.36 111.212 50.111 28.2374 19.36 12.9619])';
% 2 level oriented + packlift MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A6 = psnr([229.581 133.398 64.0705 37.3272 24.9564 17.2263])';
% HPF prelift MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A7 = psnr([220.895 129.181 65.2676 37.2316 24.6778 16.6229])';
% 2 level HPF prelift MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B) (antialiased once)
A8 = psnr([222.086 128.199 64.6631 37.6335 24.6924 16.9874])';

% Thick vertical
% Ordinary dwt MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A1 = psnr([270.906 111.382 26.4408 11.3685 3.40265 1.57077])';
% Oriented wavelet MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A3 = psnr([257.306 82.7419 17.1686 4.95542 2.37916 1.07448])';
% 2 level oriented MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A4 = psnr([233.493 68.2856 13.2213 4.29295 2.04008 0.940356])';
% HPF prelift MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A7 = psnr([276.893 101.199 21.5319 8.19918 3.38796 1.79382])';
% 2 level HPF prelift MSEs w5x3 Cdecomp:B(BH-H-:BVV--:-),B(H:V:B)
A8 = psnr([247.537 86.3365 17.6104 7.25431 2.64812 1.54708])';


figure, plot(rate,[A3 A4 A7 A8])
legend('orient','orient2','hpf prelift','hpf prelift 2 oriented levels')

legend('orient','orient2','packlift 2 oriented levels','hpf prelift','hpf prelift 2 oriented levels')
% legend('base','packlift','orient packet','orient2','orient1','orient2 + packlift')