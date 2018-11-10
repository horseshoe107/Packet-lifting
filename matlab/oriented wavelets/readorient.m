function [V,H] = readorient(fname)
% [V,H] = readorient(fname)
%
% Reads a orientation .dat file and returns two matrices:
% V is the horizontal shifts corresponding to the vertical transform
% H is the vertical shifts corresponding to the horizontal transform
fid = fopen(fname);
tline = fgetl(fid);
A = fscanf(fid,'%d%*c');
fclose(fid);

header = sscanf(tline,'%d');
w = header(1);
h = header(2);
blksz = header(3);
% prec = header(4);
% maxshft = header(5);

B = reshape(A,w/blksz*2,h/blksz)';
V = B(:,1:2:size(B,2));
H = B(:,2:2:size(B,2));