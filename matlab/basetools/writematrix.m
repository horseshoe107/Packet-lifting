function writematrix(fid,A)
% writematrix(fid,A)
%
% write the contents of A out to file fid in row-major order with 10
% decimal precision and newlines for each new row
for m=1:size(A,1)
    for n=1:size(A,2)
        if A(m,n)<0
            fprintf(fid,'%.10f ',A(m,n));
        else
            fprintf(fid,' %10.10f ',A(m,n));
        end
    end
    fprintf(fid,'\n');
end