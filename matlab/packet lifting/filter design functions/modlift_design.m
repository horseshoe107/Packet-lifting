function [Ht,Hc] = modlift_design(passband,passband_ripple,stopband, ...
                             stopband_ripple,forced_zeros,max_cancel_mag)
% [Ht,Hc] = modlift(P_freq,P_ripple,S_freq,S_ripple,forced_zeros,max_cancel_mag)
% Generates three filters, such that Qz_max is the cancel
% filter and Pz, convolved with Qz_min is the transfer
% filter.  Pz is the unit gain filter whose zeros are at the frequencies described
% by the `forced_zeros' row vector.  These zeros always appear in pairs; normally,
% you will only place one pair of such zeros, the normalized frequency 1 (Nyquist).
% The P_freq and S_freq arguments identify the normalized frequencies (1 is Nyquist)
% of the passband and stopband.  The P_ripple and S_ripple arguments identify the
% corresponding ripple tolerances.  Finally, `max_cancel_mag' identifies the
% maximum magnitude of the cancel filter at any frequency.  This constraint may
% not be compatible with the others, in which case the function will exit with
% an error message.

% this is the function that will run
max_Qz_mag = max_cancel_mag^2;

% find the forced component (Pz contributes to Ht only)
Pz = 1;
for n=1:size(forced_zeros,2)
  theta = forced_zeros(n)*pi;
  conj_pair = [1 -2*cos(theta) 1];
  Pz = conv(Pz,conj_pair) / sum(conj_pair);
end

% Design criterion for Q(z) - find passband bounds
passband_segs = 16; % change to adjust precision
num_band_segs = passband_segs+1; % unnecessary?
freqs = zeros(1,2*num_band_segs);
upper_bound = zeros(1,2*num_band_segs);
lower_bound = zeros(1,2*num_band_segs);
for n=0:(passband_segs-1)
  freqs(1+2*n) = passband*(n+0.1)/passband_segs; % get 2 frequencies for each segment
  freqs(2+2*n) = passband*(n+0.9)/passband_segs;
  Pz_mag = abs(exp(-j*freqs(1+2*n)*pi*(1:size(Pz,2)))*Pz');
  upper_bound(1+2*n) = (1+passband_ripple)/(Pz_mag+eps);% - stopband_ripple;
  lower_bound(1+2*n) = (1-passband_ripple)/(Pz_mag+eps);% - stopband_ripple;
  Pz_mag = abs(exp(-j*freqs(2+2*n)*pi*(1:size(Pz,2)))*Pz');
  upper_bound(2+2*n) = (1+passband_ripple)/(Pz_mag+eps);% - stopband_ripple;
  lower_bound(2+2*n) = (1-passband_ripple)/(Pz_mag+eps);% - stopband_ripple;
end
upper_bound = min(upper_bound,max_Qz_mag);
for n=1:(2*passband_segs)
  if (upper_bound(n) < lower_bound(n)+0.001)
    error('Upper bound on the cancel filter magnitude is too tight');
  end
end

% set stopband bounds
freqs(2*num_band_segs-1:2*num_band_segs) = [stopband 1.0];
upper_bound(2*num_band_segs-1:2*num_band_segs) = stopband_ripple;
lower_bound(2*num_band_segs-1:2*num_band_segs) = -stopband_ripple;

% establish target and weight functions from upper and lower bounds
target = (upper_bound + lower_bound)/2;
Delta = upper_bound - lower_bound;
Delta1 = Delta(1:2:(2*num_band_segs));
Delta2 = Delta(2:2:(2*num_band_segs));
weights = 1./min(Delta1,Delta2);

% pre-adjust target
target(1:(2*passband_segs)) = target(1:(2*passband_segs)) - stopband_ripple/2;
% produces final passband gain ~< 1 since true ripple adjustment (below)
% will always be less than stopband_ripple/2

% Use method of parks-mclellan to design Qz
found = 0;
N = 2;
while ((found == 0) & (N < 60))
  N = N+2;
  Qz = firpm(N,freqs,target,weights);
  found = 1;
  w=(1:1000)*pi/1000; % change to adjust precision of estimate
  tmp_freqs = real(exp(-j*w'*((-N/2):(N/2)))*Qz');
  min_val = min(tmp_freqs);
  Qz(1+N/2) = Qz(1+N/2) - min_val;
  for n=1:(2*num_band_segs)
    Qz_freqs(n) = real(exp(-j*freqs(n)*pi*((-N/2):(N/2)))*Qz');
    if (Qz_freqs(n) > upper_bound(n))
      found = 0;
    elseif (Qz_freqs(n) < lower_bound(n))
      found = 0;
    end
  end
end
disp(['N = ' num2str(N)])

% feedback - redesign
target(1:(2*passband_segs)) = target(1:(2*passband_segs)) + stopband_ripple/2 + min_val;
Qz = firpm(N,freqs,target,weights);
Qz(1+N/2) = Qz(1+N/2) - min_val;

% factorise into min/max phase components
Qz_roots = roots(Qz);
num_roots = size(Qz_roots,1);
Qz_mag = abs(Qz_roots);
minp_indices = find(Qz_mag < 1);
maxp_indices = find(Qz_mag > 1);
% find stopband zeros and move to unit circle
stopband_rounding = 0.05;
lower_thresh = 1 - stopband_rounding;
upper_thresh = 1/lower_thresh;
for n=1:num_roots
    Rmag = abs(Qz_roots(n));
    if (lower_thresh < Rmag) && (Rmag < upper_thresh)
        Qz_roots(n) = exp(j*angle(Qz_roots(n)));
    end
end

Qz_minp = real(poly(Qz_roots(minp_indices)));
Qz_maxp = real(poly(Qz_roots(maxp_indices)));
Qz_minp = Qz_minp * sqrt(sum(Qz))/sum(Qz_minp);
Qz_maxp = Qz_maxp * sqrt(sum(Qz))/sum(Qz_maxp);

Ht = conv(Pz,Qz_minp);
Hc = Qz_maxp;