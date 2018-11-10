function Ain = rawio(fname,h,w)
% Ain = rawio(fname,h,w)
%       rawio(fname,Aout)
%
% Read/write with file fname. On writing out, Aout must have values in the
% range [-2^9,2^9-1] to avoid overflow errors. Numbers are multiplied by
% 2^6 before being rounded off and written out as 16 bit signed integers.
% Read operation is the reverse of this process.
exp = 6;
if nargin==2
    if (min(min(h))<-2^9)||(max(max(h))>=2^9)
        error('Data must be in the range [-2^9,2^9-1]')
    end
    fid = fopen(fname,'wb');
    fwrite(fid,h'*2^exp,'int16');
else
    fid = fopen(fname,'rb');
    Ain = fread(fid,[w,inf],'int16')'*2^(-exp);
    if size(Ain,1)~=h
        error('File does not match dimensions given')
    end
end
fclose(fid);