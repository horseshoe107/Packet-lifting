passband = 0.45;
passband_ripple = 0.3;  % Applies to the composite filter
stopband = 0.575;
stopband_ripple = 0.1; % Applies to the composite filter
max_cancel_mag = 1.3;

% Specify forced zeros for the combined filter via their
% normalized frequencies (0 to 1).  These zeros will occur in
% conjugate pairs, regardless of their frequencies.  The remainder
% of the filter design problem consists of finding a linear phase
% common part, Q(z), which can be factorized into minimum and
% maximum phase filters.  We do not include the effect of the
% `forced_zeros' on the stopband, specifically so that we can get
% a stopband whose zeros all exactly touch the origin.  It also
% means that the `forced_zeros' provide stronger attenuation in
% the stopband.
forced_zeros = [1]; % Can be a row vector with more entries

[Ht,Hc] = modlift_design(passband,passband_ripple, ...
                                      stopband,stopband_ripple,...
                                      forced_zeros,max_cancel_mag);

% filttest1
[Ht,Hc] = modlift_design(0.45,0.1,0.65,0.1,[1],50);
% filttest3
[Ht,Hc] = modlift_design(0.4,0.1,0.6,0.05,[1],50);
Hp = packlift_preupdate(Ht,Hc);
% filttest4 - low order
[Ht,Hc] = modlift_design(0.4,0.1,0.6,0.1,[1],50);

% w = (0:100)'*pi/100;
% Pz_freq = exp(-j*w*(1:size(Pz,2)))*Pz';
% Tz_freq = exp(-j*w*(1:size(Tz,2)))*Tz';
% Qz_minp_freq = exp(-j*w*(1:size(Qz_minp,2)))*Qz_minp';
% Qz_maxp_freq = exp(-j*w*(1:size(Qz_maxp,2)))*Qz_minp';
% Ht_freq = exp(-j*w*(1:size(Ht,2)))*Ht';
% Hc_freq = exp(-j*w*(1:size(Hc,2)))*Hc';
% Ht_eff_freq = exp(-j*w*(1:size(Ht_eff,2)))*Ht_eff';
% figure
% plot(w/pi,abs(Ht_freq),'b'); hold on;
% plot(w/pi,abs(Hc_freq),'r');
% plot(w/pi,abs(Tz_freq),'m');
% plot(w/pi,abs(Ht_eff_freq),'b--');
% plot([passband,passband+eps],[0,1],'k--');
% plot([stopband,stopband+eps],[0,1],'k--');
% plot([0.5,0.5+eps],[0,1],'g--'); hold off;
% legend('ht','hc','tz','ht_eff')