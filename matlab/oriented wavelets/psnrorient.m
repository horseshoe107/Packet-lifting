R = [0.1 0.2 0.4 0.6 0.8 1];

% BARBARA - Cdecomp=B(BH-H-:BVV--:-),B(B:B:B),B(-:-:-), 5x3 dwt base
A1 = psnr([238.491 138.475 67.1136 37.4123 24.7742 16.8553]); % baseline
A2 = psnr([247.094 146.898 77.1338 46.2092 31.6937 22.078]);  % antialised
A3 = psnr([198.533 106.331 47.6029 27.7785 18.8514 12.828]);  % oriented (barb4.dat)
A4 = psnr([201.763 111.035 54.5271 33.5142 24.0046 16.1946]); % aa+orient (barb4)
A5 = psnr([209.157 115.688 53.5402 29.4871 20.061 13.4033]);  % 2*orient (barb4) - blksize halved

% BARBARA - Cdecomp=B(BH-H-:BVV--:-),B(B:B:B),B(-:-:-), 9x7 dwt base
B1 = psnr([194.415 105.855 47.6255 27.3545 17.4965 11.5959]); % baseline
B2 = psnr([200.555 114.944 56.2522 33.3071 21.8298 14.9103]); % antialiased
B3 = psnr([173.127 88.4925 39.244 22.4439 15.0302 10.6979]);  % oriented (barb4.dat)
B4 = psnr([180.981 95.257 48.0732 27.4313 18.8727 13.426]);   % aa+orient (barb4)
B5 = psnr([176.239 91.2995 40.5399 22.8648 15.249 10.8316]);  % 2*orient (barb4) - option 2

figure(1), plot(R,[A1' A2' A3' A4']), grid on
% legend('Non-oriented DWT','Oriented (min energy)','Oriented (improved)',...
%     'Non-oriented antialiased','Oriented antialiased','Location','SouthEast')
figure(2), plot(R,[B1' B2' B3' B4' B5']), grid on

B6 = psnr([205.228 116.995 54.5427 31.5864 21.2592 13.9131]); % adaptive antialiasing
B6 = psnr([203.944 114.586 56.3208 31.5489 21.177 13.8284]); % reference adaptive
B7 = psnr([202.872 113.467 55.7701 31.3288 21.006 13.7234]); % perfect knowledge

% no transfer, only cancellation (barbara)
B8 = psnr([200.282 108.435 49.714 28.4574 18.2957 12.4498]);