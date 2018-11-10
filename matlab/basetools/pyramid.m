function [O1,O2] = pyramid(mode,varargin)
% [c,d] = pyramid('flierl',x)
% x = pyramid('flierl',c,d)
% [c,d] = pyramid('tran',x)
% x = pyramid('tran',c,d)
%
% Perform one level of laplacian pyramid analysis or synthesis. The mode
% argument selects:
% 'F' - 3x2 lifted pyramid as proposed by Flierl et al
% 'T' - 2x3 lifted pyramid as proposed by Tran et al
% 
% Analysis or synthesis is inferred from the number of arguments.
if any(upper(mode)=='F')
    Asteps=3;
    Ssteps=2;
else
    Asteps=2;
    if upper(mode)=='T'
        Ssteps=3;
    else % normal laplacian pyramid
        Ssteps=1;
    end
end

mpegB = [2,0,-4,-3,5,19,26,19,5,-3,-4,0,2]'/64; % passband gain 1
% h264 = [1 -5 20 20 -5 1]/32 * 2; % h264 half-pel interpolator
% modified kernel avoids half-pixel shift
h264 = [1,0,-5,0,20,32,20,0,-5,0,1]'/64 * 2; % passband gain 2

if nargin == 2 % analysis
    x = varargin{1};
    [h,w]=size(x);
    % filter with mpegB and downsample
    a = filtextend(filtextend(x,mpegB,1)',mpegB,1)';
    c = a(1:2:h,1:2:w); % prelimary coarse band
    % upsample, filter with h264 and subtract
    b = zeros(h,w); b(1:2:h,1:2:w) = c;
    b = filtextend(filtextend(b,h264,1)',h264,1)';
    d = x - b; % difference band
    if Asteps==3 % update step: filter with mpegB, downsample and add
        e = filtextend(filtextend(d,mpegB,1)',mpegB,1)';
        c = c + e(1:2:h,1:2:w);
    end
    % cleanup
    O1 = c; O2 = d;
end

if nargin == 3 % synthesis
    c = varargin{1};
    d = varargin{2};
    [h,w] = size(d);
    if Ssteps==3 % upsample, filter with h264 and add
        b = zeros(h,w); b(1:2:h,1:2:w) = c;
        b = filtextend(filtextend(b,h264,1)',h264,1)';
        d = d + b;
    end
    if Ssteps>=2 % filter with mpegB, downsample and subtract
        a = filtextend(filtextend(d,mpegB,1)',mpegB,1)';
        c = c - a(1:2:h,1:2:w);
    end
    % upsample, filter with h264 and add
    b = zeros(h,w); b(1:2:h,1:2:w) = c;
    b = filtextend(filtextend(b,h264,1)',h264,1)';
    O1 = d + b;
end