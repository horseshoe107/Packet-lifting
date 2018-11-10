% plot the frequency response for pairs of analysis filters
% w = 0:0.01:1;
% h0 = h0_97; h1=h1_97;
% a = h0*exp(-j*pi*(1:length(h0))'*w);
% b = h1*exp(-j*pi*(1:length(h1))'*w);
% figure(1), plot(w,abs(a),w,abs(b),'LineWidth',2)
% axis([0 1 0 1.2]), grid on
% ylabel('Magnitude'), xlabel('Normalised Frequency (\times \pi)')
% legend('Low pass\newlineanalysis vector',...
% 'High pass\newlineanalysis vector','Location','West')
% set(gca,'Xtick',[0 0.25 0.5 0.75 1])

Y = double(imread('city0.pgm'));
[h,w] = size(Y);
Ht = [0.01844836   0.00352212  -0.05923167  -0.02809127   0.05832551...
    0.01403725  -0.11084035  -0.06290821  0.21398923  0.43670698  0.35060415  0.10802837];
Hc = [0.43230247   0.53774621   0.23867446  -0.15904231  -0.17211161...
    0.06003198   0.10818444  -0.04334815  -0.13367273   0.07382572];
Hp = [0.35528428   0.47972674   0.23130024  -0.07246351  -0.16015160...
    -0.02395411   0.08505587   0.01536558  -0.03851795  -0.11773997];
mode = '';
Tgrid = ones(h/4,w/4);

X = analysis_9_7(analysis_9_7(Y',1,0)',1,0);
X = analysis_9_7(X,2,0);
X = analysis_9_7(X,2,1);
X = analysis_9_7(X',2,0)';
X = analysis_9_7(X',2,1)';

% define the swap bands
LL_HL = X(1:4:h,3:4:w);
LL_LH = X(3:4:h,1:4:w);
LL_HH = X(3:4:h,3:4:w);
HL_LL = X(1:4:h,2:4:w);
HL_LH = X(3:4:h,2:4:w);
LH_LL = X(2:4:h,1:4:w);
LH_HL = X(2:4:h,3:4:w);

% vertical first
[LL_LHpost,LH_LLpost] = packswap(LL_LH,LH_LL,Ht,Tgrid,[mode 'A'],fliplr(Hc),fliplr(Hp));

% file contains energies in the transfer bands on analysis and also for
% synthesis for all 6 compression rates
E = csvread('cancel_energy.dat'); W=size(E,2);
E = csvread('energiesfull.dat'); W=size(E,2);
for n=1:4
    if n<3  h=115; w=126;
    else    h=126; w=115;
    end
    htE = reshape(E(n,1:2:W)',w,h)';
    donorE = reshape(E(n,2:2:W)',w,h)';
    if n>2 n=n-4; end % synthesis reverses vertical packets with horizontal
    htE_10bpp = reshape(E(n+26,1:2:W)',w,h)';
    donorE_10bpp = reshape(E(n+26,2:2:W)',w,h)';
    htE_08bpp = reshape(E(n+22,1:2:W)',w,h)';
    donorE_08bpp = reshape(E(n+22,2:2:W)',w,h)';
    htE_06bpp = reshape(E(n+18,1:2:W)',w,h)';
    donorE_06bpp = reshape(E(n+18,2:2:W)',w,h)';
    htE_04bpp = reshape(E(n+14,1:2:W)',w,h)';
    donorE_04bpp = reshape(E(n+14,2:2:W)',w,h)';
    htE_02bpp = reshape(E(n+10,1:2:W)',w,h)';
    donorE_02bpp = reshape(E(n+10,2:2:W)',w,h)';
    htE_01bpp = reshape(E(n+6,1:2:W)',w,h)';
    donorE_01bpp = reshape(E(n+6,2:2:W)',w,h)';
    F(1)=mean(mean(htE));
    F(2)=mean(mean(htE_10bpp));
    F(3)=mean(mean(htE_08bpp));
    F(4)=mean(mean(htE_06bpp));
    F(5)=mean(mean(htE_04bpp));
    F(6)=mean(mean(htE_02bpp));
    F(7)=mean(mean(htE_01bpp));
    G(1)=mean(mean(donorE));
    G(2)=mean(mean(donorE_10bpp));
    G(3)=mean(mean(donorE_08bpp));
    G(4)=mean(mean(donorE_06bpp));
    G(5)=mean(mean(donorE_04bpp));
    G(6)=mean(mean(donorE_02bpp));
    G(7)=mean(mean(donorE_01bpp));
    disp(F), disp(G), disp(F./G)
end

figure, imshow(htE./(donorE+0.1),[0 2])
figure, imshow(htE_10bpp./(donorE_10bpp+0.1),[0 2])
figure, imshow(htE_08bpp./(donorE_08bpp+0.1),[0 2])
figure, imshow(htE_04bpp./(donorE_04bpp+0.1),[0 2])
figure, imshow(htE_02bpp./(donorE_02bpp+0.1),[0 2])

figure, imshow(htE./(donorE+1000),[0 2])
figure, imshow(htE_10bpp./(donorE_10bpp+1000),[0 2])
figure, imshow(htE_08bpp./(donorE_08bpp+1000),[0 2])
figure, imshow(htE_04bpp./(donorE_04bpp+1000),[0 2])
figure, imshow(htE_02bpp./(donorE_02bpp+1000),[0 2])

A = csvread('alpha_cancel.dat');
csvwrite('alpha_cancel.dat',A(:,1:14490));
h=106; w=126; n=2; % cancellation step
h=115; w=126; n=2; % transfer step
A1 = reshape(A(n,:)',w,h)';
if n>2 n=n-4; end
A1_01bpp = reshape(A(n+6,:)',w,h)';
A1_02bpp = reshape(A(n+10,:)',w,h)';
A1_04bpp = reshape(A(n+14,:)',w,h)';
A1_06bpp = reshape(A(n+18,:)',w,h)';
A1_08bpp = reshape(A(n+22,:)',w,h)';
A1_10bpp = reshape(A(n+26,:)',w,h)';
figure, imshow(A1,[])
figure, imshow(A1_10bpp,[])
figure, imshow(A1_08bpp,[])
figure, imshow(A1_04bpp,[])
figure, imshow(A1_02bpp,[])