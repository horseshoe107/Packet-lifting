function out = removepadding(M)
% out = removepadding(M)
%
% remove zero padding from matrix M
% usage example:
% removepadding([0 0 0 0 0;
%                0 0 0 0 0;
%                0 1 4 6 0;
%                0 0 5 3 0;
%                0 0 0 0 0])
% ans = 1 4 6
%       0 5 3
[h,w] = size(M);
T = (M==0);
Ystart=h;
for y=1:(h-1)
    if all(T(y,:))
        continue
    else
        Ystart=y;
        break;
    end
end
Yend=1;
for y=h:-1:(Ystart+1)
    if all(T(y,:))
        continue
    else
        Yend=y;
        break
    end
end
Xstart=w;
for x=1:(w-1)
    if all(T(:,x))
        continue
    else
        Xstart=x;
        break;
    end
end
Xend=1;
for x=w:-1:(Xstart+1)
    if all(T(:,x))
        continue
    else
        Xend=x;
        break
    end
end
out = M(Ystart:Yend,Xstart:Xend);