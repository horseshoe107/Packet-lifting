function Y = interleave(x1,x2)
% interleave 2 sequences together. x1 should have an odd number of samples,
% while x2 should have an even number of samples

len1 = length(x1);
len2 = length(x2);
diff = abs(len1-len2);
if mod(len1,2)==mod(len2,2) % error checking
    error(['The two sequences cannot both have odd/even samples'...
        ' - interleaving is undefined']);
end
if len1 > len2
    Y(1:2:(2*len1-1))=x1;
    Y((1+diff):2:(2*len2-1+diff))=x2;
else
    Y(1:2:(2*len2-1))=x2;
    Y((1+diff):2:(2*len1-1+diff))=x1;
end