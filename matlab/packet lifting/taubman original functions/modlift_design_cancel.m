function Hc_opt = modlift_design_cancel(Ht,Hc,Hp,alpha)
% Hc_opt = modlift_design_update(Ht,Hc,Hp)
% Least squares filter design procedure to find an improved
% cancellation filter, Hc, assuming that the transfer (Ht)
% and pre-update (Hp) filters are fixed.  To understand the
% formulation, we consider two ways in which Hc has an impact
% on the performance of our modulated lifting system.
%
% We begin by consider quantization noises, N1 and N2, during
% synthesis, corresponding to the donor (low-pass) and
% receiver (high-pass) packet subbands, respectively.  After
% application of the synthesis system, these noise powers
% become:
%    N1' = N1 * (I - HpHt) + N2 * A
% and
%    N2' = N1 * (-Ht) + N2 * B
% where the filters are considered as composed operators,
% I is the identity operator, and * denotes application of
% the filtering operator to a noise signal.  The A and B
% filters contain the dependence on the Hc filter which we
% wish to design.  Specifically,
%    A = Hc + Hp - HpHtHc
% and
%    B = I - HtHc
% From the perspective of the synthesis system, we can say
% the following.  If the synthesis basis functions end up being
% close to orthonormal (a desirable situation, which our
% optimization tends to encourage), noise signals N1 and N2
% should be selected to have roughly equal powers by a rate
% distortion optimization framework, at least at very high
% bit-rates.  We also assume that the noise processes are
% both white (again, this is reasonable at very high bit-rates).
% In this case, our optimization objective is to minimize
% the total power gain associated with operators A and B.
% Equivalently, we wish to minimize
%    Js = \int | \hat{A}(\omega) |^2 + | \hat{B}(\omega) |^2 d\omega
%
% We turn our attention now to the impact on analysis, considering
% donor and receiver subbands whose content is initially described
% by signals E1 and E2, respectively.  After application of the
% analysis system, these signals become
%    E1' = E1 * B  -  E2 * A
% and
%    E2' = E1 * Ht + E2 * (I - HtHp)
% where A and B are the same operators described above.  Since
% the objective of the cancellation step is to minimize the
% energy of signal E1', it is again reasonable to minimize the
% total power gain associated with operators A and B.  We should
% consider the power spectra of E1 and E2, which are neither equal
% nor uniform -- unlike the quantization noise spectra considered
% above.  To this end, we use a simple model of the power spectrum
% of natural images, of the form \Gamma(\omega) = K / |\omega|^2.
% Since the subbands which are being operated on are created by
% a wavelet packet transform, we need to be careful in correctly
% identifying the value of \omega here.  Specifically, the spectrum
% of E1, in its packet subband can be identified as
%    \Gamma_E1(\omega) = K / (\pi - |\omega|/2)^2
% while that of E2 may be identified as
%    \Gamma_E2(\omega) = K / (2 \pi - |\omega|/2)^2
% If we focus only on the minimization of the energy E1', our
% optimization objective would be
%    Ja = \int | \hat{A}(\omega) |^2 / (2 \pi - |\omega|/2)^2 +
%              | \hat{B}(\omega) |^2 / (\pi - |\omega|/2)^2   d\omega
%
% Since we do not know the value of K, or even how to decide
% between the relative importance of objectives Ja and Js, we
% opt simply to minimize a weighted sum of these to objectives.
% Specifically, we find the Hc filter which minimizes
%    Js + 2 * pi * alpha * Ja.
% Smaller values of alpha tend to ignore the differences in
% power spectra associated with E1 and E2, whereas larger values
% of alpha tend to emphasize the cancellation of original signal
% content, E1, in the donor band.
%
% To recap, we choose Hc which minimizes the frequency weighted
% power gains associated with the two filters:
%    A = Hc + Hp - HpHtHc  =  Hp - Hc(HtHp - I)
% and
%    B = I - HtHc
%
% If the `alpha' argument is negative, the present function has the
% special behaviour of just returning the supplied Hc value as-is,
% without further optimization

Hc_opt = Hc
if (alpha < 0.0)
   return;
end

Nc = -floor((length(Hc)-1)/2);
Pc = floor(length(Hc)/2);
Nt = -floor(length(Ht)/2);
Pt = floor((length(Ht)-1)/2);
if ((Nc+Nt+Pc+Pt) ~= 0)
   error('Ht and Hc filters have orders which are incompatible with end-to-end zero phase');
end
Np = -floor((length(Hp)-1)/2);
Pp = floor(length(Hp)/2);

w = (0:2047)*pi/1024;
abs_w = w;
abs_w(1025:2048) = abs_w(1024:-1:1);
Hc_freq = exp(-j*w'*(Nc:Pc))*Hc';
Ht_freq = exp(-j*w'*(Nt:Pt))*Ht';
Hp_freq = exp(-j*w'*(Np:Pp))*Hp';
Res_freq = Hp_freq .* Ht_freq - 1;  % (HtHp - I)
E1_freq = 1.0 ./ (pi - abs_w'/2).^2;
E2_freq = 1.0 ./ (2*pi - abs_w'/2).^2;

alpha = alpha * 2 * pi;

desired_freq = conj(Ht_freq) .* (1+alpha*E1_freq) ...
             + Hp_freq .* conj(Res_freq) .* (1+alpha*E2_freq);
weight_freq = Ht_freq .* conj(Ht_freq) .* (1+alpha*E1_freq) ...
             + Res_freq .* conj(Res_freq) .* (1+alpha*E2_freq);
R_cor = zeros(1,(Pc-Nc+1));
D_vec = zeros(Pc-Nc+1,1);
for n=1:(Pc-Nc+1)
   d_pos = n-1+Nc;
   D_vec(n) = real(exp(j*w*d_pos)*desired_freq) / length(w);   
   r_pos = n-1;
   R_cor(n) = real(exp(j*w*r_pos)*weight_freq) / length(w);
end

R_mat = zeros((Pc-Nc+1),(Pc-Nc+1));
for n=1:(Pc-Nc+1)
   for m=n:(Pc-Nc+1)
      R_mat(m,n) = R_cor(m-n+1);
      R_mat(n,m) = R_cor(m-n+1);
   end
end
Hc_opt = (inv(R_mat)*D_vec)';
Hc_opt_freq = exp(-j*w'*(Nc:Pc))*Hc_opt';


% Hp - Hc(HtHp - I)
% B = I - HtHc
%    N1' = N1 * (I - HpHt) + N2 * A
% and
%    N2' = N1 * (-Ht) + N2 * B

A_freq = Hp_freq - Hc_opt_freq .* Res_freq;
B_freq = 1 - Ht_freq .* Hc_opt_freq;
A_energy = sum(abs(A_freq).^2) / length(w);
B_energy = sum(abs(B_freq).^2) / length(w);
N1_energy = A_energy + sum(abs(Res_freq).^2)/length(w);
N2_energy = B_energy + sum(abs(Ht_freq).^2)/length(w);
J = sum(abs(A_freq).^2 .* (1+alpha*E2_freq) ...
        + abs(B_freq).^2 .* (1+alpha*E1_freq)) / length(w);

Aorig_freq = Hp_freq - Hc_freq .* Res_freq;
Borig_freq = 1 - Ht_freq .* Hc_freq;
Aorig_energy = sum(abs(Aorig_freq).^2) / length(w);
Borig_energy = sum(abs(Borig_freq).^2) / length(w);
N1orig_energy = Aorig_energy + sum(abs(Res_freq).^2)/length(w);
N2orig_energy = Borig_energy + sum(abs(Ht_freq).^2)/length(w);
Jorig = sum(abs(Aorig_freq).^2 .* (1+alpha*E2_freq) ...
            + abs(Borig_freq).^2 .* (1+alpha*E1_freq)) / length(w);

% Write comparison statistics
Jvals = [Jorig, J]
AB_energies = [Aorig_energy  A_energy;
               Borig_energy  B_energy;
               Aorig_energy+Borig_energy  A_energy+B_energy]
N1N2_energies = [N1orig_energy  N1_energy;
                 N2orig_energy  N2_energy;
                 N1orig_energy+N2orig_energy  N1_energy+N2_energy]
              

