function [Ht_dm,Hc_dm,Hp_dm, Ht_em,Hc_em,Hp_em] = ...
   modlift_generate(passband,passband_ripple, ...
                    stopband,stopband_ripple, ...
                    forced_zeros,max_cancel_mag,alpha)
% [Ht_dm,Hc_dm,Hp_dm,Ht_em,Hc_em,Hp_em] =
%     modlift_generate(P_freq,P_ripple,S_freq,S_ripple,forced_zeros,max_cancel_mag,alpha)
% Generates the transfer, cancellation and pre-update filters, (Ht_dm, Hc_dm and Hp_dm),
% as well as equal-magnitude versions of these filters (Ht_em, Hc_em, Hp_em).  Also plots
% their responses and writes filter files named "filters_dm.dat" and "filters_em.dat"
% which can be used directly by the "modlift.exe" program.  The `alpha' argument
% is passed into `modlift_design_cancel'.  It controls the relative significance given
% to the cancellation of original content in the donor (lower frequency) subband, as
% opposed to the minimization of aliasing energy introduced from the receiver (higher
% frequency) subband.  For more on this parameter, consult the description of
% `modlift_design_cancel'.

[Qz_minp,Qz_maxp,Pz] = modlift_design(passband,passband_ripple, ...
                                      stopband,stopband_ripple,...
                                      forced_zeros,max_cancel_mag);
Ht = conv(Pz,Qz_minp);
Hc = Qz_maxp;
Tz = conv(Ht,Hc);
Tz_delta = zeros(1,size(Tz,2));
Tz_delta((size(Tz,2)+1)/2) = 1;
Ht_eff = conv((2*Tz_delta - Tz),Ht);
Hc_eff = conv((2*Tz_delta - Tz),Hc);

Hp = modlift_design_preupdate(Ht,Hc);
Hp_eff = conv((Tz_delta - Tz),Hp) + conv(Tz_delta,Hc);
Hc_opt = modlift_design_cancel(Ht,Hc,Hp,alpha);

w = (0:100)'*pi/100;
Tz_freq = exp(-j*w*(1:size(Tz,2)))*Tz';
Ht_freq = exp(-j*w*(1:size(Ht,2)))*Ht';
Hc_freq = exp(-j*w*(1:size(Hc,2)))*Hc';
Hc_opt_freq = exp(-j*w*(1:size(Hc_opt,2)))*Hc_opt';
Ht_eff_freq = exp(-j*w*(1:size(Ht_eff,2)))*Ht_eff';
Hc_eff_freq = exp(-j*w*(1:size(Hc_eff,2)))*Hc_eff';
Hp_eff_freq = exp(-j*w*(1:size(Hp_eff,2)))*Hp_eff';

figure;
plot(w/pi,abs(Ht_freq),'b'); hold on;
plot(w/pi,abs(Hc_freq),'r');
plot(w/pi,abs(Hc_opt_freq),'g');
plot(w/pi,abs(Tz_freq),'m');
plot(w/pi,abs(Ht_eff_freq),'b--');
plot(w/pi,abs(Hc_eff_freq),'r--');
plot(w/pi,abs(Hp_eff_freq),'c');
plot([passband,passband+eps],[0,1],'k--');
plot([stopband,stopband+eps],[0,1],'k--');
plot([0.5,0.5+eps],[0,1],'g--'); hold off;

Hc_len = length(Hc);
Ht_len = length(Ht);
Hp_len = length(Hp);
fid = fopen('filters_dm.dat','w');
fprintf(fid,'%d',Ht_len);
for n=1:Ht_len
  fprintf(fid,' %12.8f',Ht(n));
end
fprintf(fid,'\n');
fprintf(fid,' %d',Hc_len);
for n=1:Hc_len
  fprintf(fid,' %12.8f',Hc_opt(n));
end
fprintf(fid,'\n');
fprintf(fid,' %d',Hp_len);
for n=1:Hp_len
  fprintf(fid,' %12.8f',Hp(n));
end
fprintf(fid,'\n');
fclose(fid);

Hc_dm = Hc;
Ht_dm = Ht;
Hp_dm = Hp;
Ht_em = [0];
Hc_em = [0];
Hp_em = [0];

if ((length(forced_zeros) ~= 1) | (forced_zeros(1) ~= 1.0))
   return;
end

Ht = conv([0.5,0.5],Qz_minp);
Hc = conv([0.5,0.5],Qz_maxp);
Tz = conv(Ht,Hc);
Tz_delta = zeros(1,size(Tz,2));
Tz_delta((size(Tz,2)+1)/2) = 1;
Hc_eff = conv((2*Tz_delta - Tz),Hc);

Hp = modlift_design_preupdate(Ht,Hc);
Hp_eff = conv((Tz_delta - Tz),Hp) + conv(Tz_delta,Hc);

Tz_freq = exp(-j*w*(1:size(Tz,2)))*Tz';
Ht_freq = exp(-j*w*(1:size(Ht,2)))*Ht';
Hc_freq = exp(-j*w*(1:size(Hc,2)))*Hc';
Hc_eff_freq = exp(-j*w*(1:size(Hc_eff,2)))*Hc_eff';
Hp_eff_freq = exp(-j*w*(1:size(Hp_eff,2)))*Hp_eff';

figure;
plot(w/pi,abs(Ht_freq),'b'); hold on;
plot(w/pi,abs(Hc_freq),'r');
plot(w/pi,abs(Tz_freq),'m');
plot(w/pi,abs(Hc_eff_freq),'r--');
plot(w/pi,abs(Hp_eff_freq),'c');
plot([passband,passband+eps],[0,1],'k--');
plot([stopband,stopband+eps],[0,1],'k--');
plot([0.5,0.5+eps],[0,1],'g--'); hold off;

Hc_len = length(Hc);
Ht_len = length(Ht);
Hp_len = length(Hp);
fid = fopen('filters_em.dat','w');
fprintf(fid,'%d',Ht_len);
for n=1:Ht_len
  fprintf(fid,' %12.8f',Ht(n));
end
fprintf(fid,'\n');
fprintf(fid,' %d',Hc_len);
for n=1:Hc_len
  fprintf(fid,' %12.8f',Hc(n));
end
fprintf(fid,'\n');
fprintf(fid,' %d',Hp_len);
for n=1:Hp_len
  fprintf(fid,' %12.8f',Hp(n));
end
fprintf(fid,'\n');
fclose(fid);

Hc_em = Hc;
Ht_em = Ht;
Hp_em = Hp;