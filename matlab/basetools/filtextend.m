function Y = filtextend(X,Fnum,Fden,bext,boundary)
% Y = filtextend(X,Fnum,Fden,bext*,boundary*)
%
% Filter the columns of X with coefficients defined by Fnum and Fden. Note
% for odd order FIR filters, the signal will be shifted forward by a
% half-pixel.
%
% bext sets the method for both FIR filtering boundary extension:
% 0-> zero extension
% 1-> symmetric extension (default mode)
% -1-> produce outputs only where full supported, placing zeros elsewhere
% and for IIR filtering, initial conditions which assume the neighbouring
% region:
% 0-> zero extension
% 1-> constant extension of bordering value (default mode)
% -1-> alternating extension of bordering value
%
% boundary sets the region over which to calculate the output. It defaults
% to [1,size(X,1)] - ie, Y will be the same size as X.
% User set boundaries cannot exceed one symmetric extension of the image X.

% flip column vectors to row vectors
if size(Fnum,2)==1, Fnum=Fnum'; end
if size(Fden,2)==1, Fden=Fden'; end
% warning disabled
% if mod(size(Fnum,2),2)==0, warning('Odd order filter; half pel shift'); end
[h,w] = size(X);
if nargin<4 % boundary extension mode
    bext = 1; end
if nargin<5 % desired output location/dimensions
    boundary = [1 h]; end
a = boundary(1); b = boundary(2);
if length(boundary)>=4
    c = boundary(3); d = boundary(4);
else c = 0; d = 0;
end

N = size(Fnum,2)-1; % filter order
N1 = ceil(N/2)+1-a; % when N is odd, overpad N1 (1/2pel shift to right)
N2 = floor(N/2)+b-h;

% shorten for negative padding
X = X(1-min(N1,0):h+min(N2,0),:);

if size(Fden,2)==1 % FIR filter
    switch bext % boundary extension
        case 0, X = [zeros(N1,w) ; X ; zeros(N2,w)];
            Y = filter(Fnum,1,X);
            Y = Y((N+1):(N+1+b-a),:);
        case -1, Y = filter(Fnum,1,X); % filter, then overwrite boundaries with 0
            Y = [zeros(N1+c,w) ; Y((N+1+c):(size(Y,1)-d),:) ; zeros(N2+d,w)];
        otherwise, X = [flipud(X(2:N1+1,:)) ; X ; flipud(X(h-N2:h-1,:))];
            Y = filter(Fnum,1,X);
            Y = Y((N+1):(N+1+b-a),:);
    end
else % IIR filter
    % still need to sort out RHS padding
    error('Padding not yet resolved');
    switch padding % initial conditions.
        case 0, % Preceding input is zero
            Zf = filtic(Fnum,Fden,zeros(1,length(Fden)),zeros(1,length(Fnum)));
        case 1, % Preceding input has constant value equal to X(1,:)
            dcg = sum(Fnum)/sum(Fden);
            Zf = filtic(Fnum,Fden,dcg*ones(1,length(Fden)),dcg*ones(1,length(Fnum)));
        case -1, % Preceding input alternates by (-1)^n * X(1,:)
            nyqg = sum(Fnum.*(-1).^(1:length(Fnum)))/sum(Fden.*(-1).^(1:length(Fden)));
            Zf = filtic(Fnum,Fden,nyqg*(-1).^(1:length(Fden)),(-1).^(1:length(Fnum)));
    end
    Y = filter(Fnum,Fden,X,Zf'*X(1,:));
end