function [Qz_minp,Qz_maxp,Pz] = modlift_design(passband,passband_ripple, ...
   stopband,stopband_ripple,forced_zeros,max_cancel_mag)
% [Qz_min,Qz_max,Pz] = modlift(P_freq,P_ripple,S_freq,S_ripple,forced_zeros,max_cancel_mag)
% Generates three filters, such that Pz is a filter with a pair of roots at each
% of the normalized frequencies given by the `forced_zeros' row vector, and
% Qz_min and Qz_max are minimum and maximum phase filters such that
% Qz_min(z) * Qz_max(z) * Pz(z) satisfies the design constraints for a composit
% filter with the indicated passband and stopband edges (all normalized
% frequencies, with the Nyquist frequency at 1).  The P_ripple argument provides
% the ripple constraint for the passband of the composit filter.  The S_ripple
% argument identifies the stopband ripple constraint for the filter
% Qz_min(z) * Qz_max(z), so that the actual stopband of the composit filter will
% generally be smaller. Finally, `max_cancel_mag' identifies the maximum magnitude
% of the potential cancellation filter, Qz_max(z) * Pz(z), at
% any frequency.  This constraint may not be compatible with the others, in
% which case the function will exit with an error message.

max_Qz_mag = max_cancel_mag^2;

Pz = 1;
for n=1:size(forced_zeros,2)
  theta = forced_zeros(n)*pi;
  conj_pair = [1 -2*cos(theta) 1];
  Pz = conv(Pz,conj_pair) / sum(conj_pair);
end

% Find the design criterion for Q(z)
passband_segs = 16;
num_band_segs = passband_segs+1;
freqs = zeros(1,2*num_band_segs);
upper_bound = zeros(1,num_band_segs);
lower_bound = zeros(1,num_band_segs);
for n=0:(passband_segs-1)
  freqs(1+2*n) = passband*(n+0.1)/passband_segs;
  freqs(2+2*n) = passband*(n+0.9)/passband_segs;
  Pz_mag = abs(exp(-j*freqs(1+2*n)*pi*(1:size(Pz,2)))*Pz');
  upper_bound(1+2*n) = (1+passband_ripple-stopband_ripple) / (Pz_mag+eps);
  lower_bound(1+2*n) = (1-passband_ripple) / (Pz_mag+eps);
  Pz_mag = abs(exp(-j*freqs(2+2*n)*pi*(1:size(Pz,2)))*Pz');
  upper_bound(2+2*n) = (1+passband_ripple-stopband_ripple) / (Pz_mag+eps);
  lower_bound(2+2*n) = (1-passband_ripple) / (Pz_mag+eps);
end

upper_bound = min(upper_bound,max_Qz_mag);
for n=1:(2*passband_segs)
  if (upper_bound(n) < lower_bound(n)+0.001)
    error('Upper bound on the cancel filter magnitude is too tight');
  end
end

freqs(2*num_band_segs-1) = stopband;
upper_bound(2*num_band_segs-1) = stopband_ripple;
lower_bound(2*num_band_segs-1) = -stopband_ripple;
freqs(2*num_band_segs) = 1.0;
upper_bound(2*num_band_segs) = stopband_ripple;
lower_bound(2*num_band_segs) = -stopband_ripple;

weights = zeros(1,num_band_segs);
target = zeros(1,num_band_segs);
for n=0:(num_band_segs-1)
  U = upper_bound(1+2*n);
  L = lower_bound(1+2*n);
  Delta1 = U-L;
  target(1+2*n) = 0.5*(U+L);
  U = upper_bound(2+2*n);
  L = lower_bound(2+2*n);
  Delta2 = U-L;
  target(2+2*n) = 0.5*(U+L);
  Delta = min(Delta1,Delta2);
  weights(n+1) = 1/Delta;
end

found = 0;
N = 2;
while ((found == 0) & (N < 60))
  N = N+2
  Qz = remez(N,freqs,target,weights);
  found = 1;
  w=(1:1000)*pi/1000;
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

Qz_roots = roots(Qz);
num_roots = size(Qz_roots,1);
minp_indices = zeros(1,num_roots/2);
maxp_indices = zeros(1,num_roots/2);
minp_count = 0;
maxp_count = 0;
for n=1:num_roots
  if (length(find(maxp_indices == n)) == 0)  % Root is not used already
    minp_count = minp_count+1;
    minp_indices(minp_count) = n;
    inv_val = (1.0 / Qz_roots(n))';
    for m=n+1:num_roots
      if (abs(Qz_roots(m)-inv_val) < 0.01)
        maxp_count = maxp_count+1;
        maxp_indices(maxp_count) = m;
        break;
      end
    end
  end
end
Qz_minp = poly(Qz_roots(minp_indices));
Qz_maxp = poly(Qz_roots(maxp_indices));
Qz_minp = Qz_minp * sqrt(sum(Qz))/sum(Qz_minp);
Qz_maxp = Qz_maxp * sqrt(sum(Qz))/sum(Qz_maxp);
