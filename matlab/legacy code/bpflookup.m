function b = bpflookup(breakpoints,f,N)
% b = bpflookup(breakpoints,f)
% Returns the value at normalised frequency f for a notional band-pass
% filter that is simply defined as follows:
% 1. The filter has zero phase response (so b(w)==b(-w))
% 2. The magnitude response is a piecewise linear trapezoid whose
% breakpoints are set by the first argument.
%
% b = bpflookup(breakpoints,f,N)
% Returns the value at normalised frequency f for a designed filter of
% order N. The firpm design function is used to generate the desired band
% pass filter.
if (f<0)
    f=-f;
end
if (any([f breakpoints]<0))||(any([f breakpoints]>1))
    error('All frequencies must be in the normalised range [0,1]');
end
if breakpoints ~= sort(breakpoints)
    error('Breakpoints must be supplied in increasing order');
end
if nargin==2
    transitionwidth(1) = breakpoints(2) - breakpoints(1);
    transitionwidth(2) = breakpoints(4) - breakpoints(3);
    if f<=breakpoints(1)
        b=0;
    else if f>=breakpoints(4)
            b=0;
        else if f<breakpoints(2)
                b=(f-breakpoints(1))/transitionwidth(1);
            else if f<breakpoints(3)
                    b=1;
                else if f<breakpoints(4)
                        b=(breakpoints(4)-f)/transitionwidth(2);
                    else b=0;
                    end
                end
            end
        end
    end
else
    A = firpm(N,[0 breakpoints 1],[0 0 1 1 0 0]);
    A_index = -(N/2):(N/2);
    b = A*exp(-1j*A_index'*f*pi);
end