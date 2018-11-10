function pgmwrite(Y,filename)
% pgmwrite(Y,filename)
%
% save data in Y to a pgm file. each value should be in [0,255]
[h,w] = size(Y);
if (max(max(Y))>255)||(min(min(Y))<0)
    error('Data must be in the range [0,255]')
end
fid = fopen(filename,'wb');
fprintf(fid,'P5\n%d %d\n255\n',w,h);
fwrite(fid,Y','uint8');
fclose(fid);