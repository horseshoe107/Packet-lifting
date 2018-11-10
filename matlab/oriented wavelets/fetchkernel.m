function K = fetchkernel(shift)

K_LUT = [ 0    0             0             1             0             0             0
-0.0057314933  0.0213902240 -0.0798294029  0.9674625206  0.1198302579 -0.0315853847  0.0084632783 
-0.0088233996  0.0329293754 -0.1228941022  0.8796881878  0.2687588056 -0.0678352739  0.0181764069 
-0.0097770136  0.0364883116 -0.1361762329  0.7514540266  0.4325511590 -0.1018238761  0.0272836254 
-0.0090991941  0.0339586548 -0.1267354249  0.5973263672  0.5973263672 -0.1267354249  0.0339586548 ];

zshift = round(shift);
shift = (shift-zshift) * 8;
index  = abs(shift)+1;
direction = shift >= 0;

if direction
    K = K_LUT(index,:);
else
    K = fliplr(K_LUT(index,:));
end
if zshift>=0 % shift centre by zshift pixels
    K = [zeros(1,2*zshift) K]; % pad left side by 2*zshift zeros
else
    K = [K zeros(1,-2*zshift)];
end