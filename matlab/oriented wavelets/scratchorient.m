R = double(imread('barbara.pgm')); % reference image
R = double(imread('test_horz_crop.pgm'));
X = analysis_5_3(analysis_5_3(R',1,0)',1,0); % analysed

% ACSV = double(csvread('barbdwt.csv')); % csv file of shiftlift analysed image
% C = synthesis_5_3(synthesis_5_3(ACSV',1,0)',1,0); % synthesis

% Acsv1D = double(csvread('barb1D.csv')); % csv file of shiftlift 1D analysed image
% C1D = synthesis_5_3(Acsv1D,1,0); % synthesis
fid = fopen(['barbout.rawl'],'rb'); D = fread(fid,[512,512],'int16')'*2^(-6); fclose(fid);
B = synthesis_5_3(synthesis_5_3(D',1,0)',1,0);

fid = fopen(['horzoutV.rawl'],'rb'); h=256; w=512;
fid = fopen(['vertoutV.rawl'],'rb'); h=512; w=256;
A = fread(fid,[w,h],'int16')'*2^(-6);
fclose(fid);

% h=576; w=704;
h=512; w=512;
fid = fopen(['out.rawl'],'rb');
B = fread(fid,[w,h],'int16')'*2^(-6);
fclose(fid);
X = analysis_5_3(analysis_5_3(B',1,0)',1,0);
LL1 = X(1:2:h,1:2:w);
HL1 = X(1:2:h,2:2:w);
LH1 = X(2:2:h,1:2:w);
HH1 = X(2:2:h,2:2:w);
figure(1), imshow(LL1,[0 255])
figure(2), imshow(HL1,[-32,32])
figure(3), imshow(LH1,[-32,32])
figure(4), imshow(HH1,[-32,32])
figure, imshow([LL1-128 HL1*4; LH1*4 HH1*4],[-128 127])
disp(mean(mean(LL1.^2))), disp(mean(mean(HL1.^2)))
disp(mean(mean(LH1.^2))), disp(mean(mean(HH1.^2)))


% reading in orientation data for barbara
M = csvread('pens_orient.csv');
blksz=M(1,3); h=M(1,1)/blksz;
% A = reshape(M(2,:),w/blksz,h/blksz)';
% B = reshape(M(3,:),w/blksz,h/blksz)';
A = M(2:h+1,:);
B = M(2+h:2*h+1,:);
% figure, imshow([A B],[-16 16]), colormap(jet);
figure, imshow([A+B],[-16 16]), colormap(jet);
% blue and red regions suggest orientation; light blue/red means little
% shifting is required to align the transform, whereas dark blue/red is a
% large amount of shifting


R = double(imread('harbour0.pgm')); % reference image
g = [1 sqrt(2) 1;0 0 0;-1 -sqrt(2) -1];
% g = [0.0968 0.1936 0.2421 0.1936 0.0968;0 0 0 0 0;...
%     -0.0968 -0.1936 -0.2421 -0.1936 -0.0968];
P = filter2(g,R); % detect horizontal lines
Q = filter2(g',R); % vertical lines
I = sqrt(P.^2+Q.^2);
figure, imshow(I,[]);
% figure, imshow(180/pi*atan(Q./(P+eps)),[-90 90]), colormap(jet);
A = atan(Q./(P+eps)); % A in range [-pi/2,pi]
Hshift = A/pi*4*8;
Vshift = (sign(A)*pi/2 - A)/pi*4*8;
for n=1:prod(size(A))
    if (A(n)<-pi/4)||(A(n)>pi/4) % (mostly) vertical edge
        Hshift(n) = 0;
    else
        Vshift(n) = 0;
    end
end
figure, imshow(I.*sign(Hshift +Vshift),[]), colormap(jet);

Eint = csvread('edgemap.csv');
Eang = csvread('edgeangle.csv');
figure, imshow(Eint,[])
figure, imshow(Eint.*Eang,[]), colormap(jet)
figure, imshow(Eang,[]),colormap(jet)
[h,w] = size(Eint);
Sint = (Eint(1:2:h,1:2:w)+Eint(2:2:h,1:2:w)...
    +Eint(1:2:h,2:2:w)+Eint(2:2:h,2:2:w))/4;
Sang = (Eint(1:2:h,1:2:w).*Eang(1:2:h,1:2:w)...
    +Eint(2:2:h,1:2:w).*Eang(2:2:h,1:2:w)...
    +Eint(1:2:h,2:2:w).*Eang(1:2:h,2:2:w)...
    +Eint(2:2:h,2:2:w).*Eang(2:2:h,2:2:w))./(4*Sint+eps);


% sets up a striped region of predominantly horizontal orientation
h=256; ratio=2; stripewidth=20; % control settings
% ratio sets orientation angle - eg ratio=2 is equivalent to 22.5 degrees
% stripewidth sets the pixel width in the vertical direction
% We will actually construct an image 2x larger, then subsample in order to
% obtain an antialiased final result.
w=ratio*2*h+12; h=2*h+12; A = zeros(h,w);
stripewidth = 2*stripewidth;
switchperiod = stripewidth*ratio;
maximum = ratio*(h-1)+(w-1);
for y=0:h-1
    for x=0:w-1
        line_equation = ratio*y+x;
        t = mod(line_equation,2*switchperiod);
        if line_equation>maximum*3/4
            A(y+1,x+1)=1;
        end
        if line_equation<maximum/2 &&...
                (t>switchperiod)&&(t<2*switchperiod)
            A(y+1,x+1)=1;
        end
    end
end
mpegB = [2,0,-4,-3,5,19,26,19,5,-3,-4,0,2]/64;
A = filter(mpegB,1,filter(mpegB,1,A)')';
B = A(13:2:h,13:2:w); % subsample down to antialised image
B = 255*max(min(B,1),0);
figure(1), imshow(B,[]);
pgmwrite(B,'gen_orient.pgm');

A = ones(100); row=11;
A(:,1:50)=10*A(:,1:50); A(:,51:100)=245*A(:,51:100);
% imshow(A,[0 255])
X = analysis_9_7(analysis_9_7(A',1,0)',1,0);
% zero all high frequency bands (testing end-to-end filter)
Y = zeros(100); Y(1:2:100,1:2:100)=X(1:2:100,1:2:100);
Y = synthesis_9_7(synthesis_9_7(Y',1,0)',1,0);
Z = packliftanalysis(Y,Ht,Hc,Hp,'');
% figure, imshow(Z(1:2:100,1:2:100),[0 255]);
figure(5), plot(Z(row,1:2:100)); hold on
% zero the low frequency band (test cancellation filter)
Y = X; Y(1:2:100,1:2:100)=zeros(50);
Y = synthesis_9_7(synthesis_9_7(Y',1,0)',1,0);
Z = packliftanalysis(Y,Ht,Hc,Hp,'');
% figure, imshow(Z(1:2:100,1:2:100)+X(1:2:100,1:2:100),[0 255]);
figure(5), plot(Z(row,1:2:100)+X(row,1:2:100),'r');
Y = synthesis_9_7(synthesis_9_7(X',1,0)',1,0);
Z = packliftanalysis(Y,Ht,Hc,Hp,'');
figure(5), plot(Z(row,1:2:100),'m');

h=256; w=512;
A = rawio('oriented wavelets\out.rawl',h,w); B = A;
B = analysis_5_3(B,1,0);
B(1:2:h,:) = analysis_5_3(B(1:2:h,:),1,0);
B = analysis_5_3(B',1,0)';
B(:,1:2:w) = analysis_5_3(B(:,1:2:w)',1,0)';
LL2 = B(1:4:h,1:4:w);
HL2 = B(1:4:h,3:4:w);
LH2 = B(3:4:h,1:4:w);
HH2 = B(3:4:h,3:4:w);
HLL = B(1:4:h,2:2:w);
HLH = B(3:4:h,2:2:w);
LHL = B(2:2:h,1:4:w);
LHH = B(2:2:h,3:4:w);
HH1 = B(2:2:h,2:2:w);
[mse(HH1) mse(HLL) mse(HLH) mse(LHL) mse(LHH)]
[mse(HL2) mse(LH2) mse(HH2)]