function runsymbolic(sigma,varargin)
frange=(0:128)/128;
if nargin==1
    for N=1:length(frange);
        [out(N,:) T1(N,:,:) T2(N,:,:) T3(N,:,:)] = ...
            symbolic(frange(N),sigma);
    end
else
    for N=1:length(frange);
        [out(N,:) T1(N,:,:) T2(N,:,:) T3(N,:,:)] = ...
            symbolic(frange(N),sigma,varargin{1});
    end
end

% frange=(60)/128;
% [out T1 T2 T3] = symbolic(frange,sigma);

% figure, plot(frange,abs(out(:,2))), grid on
% figure, plot(frange,abs(out)), grid on
figure, plot(frange,abs(out(:,2))), grid on
% figure, plot(frange,[abs(T1(:,1,1)) abs(T1(:,1,2))]), grid on