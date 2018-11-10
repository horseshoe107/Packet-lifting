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
