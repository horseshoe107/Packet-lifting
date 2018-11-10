bitrate = 1; % selected bit rate of FULL RES image
N=150;

mov = avifile('out.avi','COMPRESSION','None','FPS',15);
cd(wkdir)
H = waitbar(0); clear F;
for n = 0:(N-1)
    A = double(imread(['city' num2str(n) tag]));

    % A = analysis_9_7(analysis_9_7(A',1,0)',1,0);
    % A=filtextend(A,mpegB,1,1,[1,size(A,1)]);
    % A=filtextend(A',mpegB,1,1,[1,size(A,2)])';
    A = wavanalysis(wavanalysis(A',1,0,h0,h1)',1,0,h0,h1);
    
    A = A(1:2:576,1:2:704);
    
    F.cdata(:,:,1) = A;
    F.cdata(:,:,2) = A;
    F.cdata(:,:,3) = A;
    F.cdata = uint8(F.cdata);
    F.colormap = [];
    mov = addframe(mov,F);
    waitbar(n/N);
end
close(H); close(mov);
cd C:\Progra~1\Matlab7\work

mov = avifile('out.avi','COMPRESSION','None','FPS',15);
cd(wkdir), Nend = 150; w=704; h=576;
H = waitbar(0); clear F;
for n = 0:(Nend-1)
    % fid = fopen(['Lift_exp\lift_' num2str(n) '_' num2str(bitrate) 'bps.rawl'],'rb');
    % fid = fopen(['DWT_exp\filt_' num2str(n) '_' num2str(bitrate) 'bps.rawl'],'rb');
    A = rawio(['LiftDump\lift_' num2str(n) '.rawl'],h,w);

    A = analysis_9_7(analysis_9_7(A',1,0)',1,0);
    A = A(1:2:576,1:2:704);

    F.cdata(:,:,1) = A+128;
    F.cdata(:,:,2) = A+128;
    F.cdata(:,:,3) = A+128;
    F.cdata = uint8(F.cdata);
    F.colormap = [];
    mov = addframe(mov,F);
    waitbar(n/Nend);
end
close(H); close(mov); clear n dims Nend A H F mov
cd C:\Progra~1\Matlab7\work