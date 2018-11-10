function Hu = modlift_design_preupdate(Ht,Hc)
% Hu = modlift_design_update(Ht,Hc)
% Least squares filter design procedure to find the pre-update filter,
% Hu, which minimizes the sum of two terms:
%     J1 = \int | \hat{A}(\omega) |^2 \Gamma_{N1}(\omega) d\omega
% and
%     J2 = \int | \hat{B}(\omega) |^2 \Gamma_{N2}(\omega) d\omega
% Here, \Gamma_{N1} models the noise power distribution due to
% quantization errors in the donor band, while \Gamma_{N2} models
% the noise power distribution due to quantization errors in the
% receiver band (i.e., the higher frequency band).  The A and B
% filters are given by:
%     A[n] = \delta[n] - Hu[n] \conv Ht[n]
% and
%     B[n] = Hc[n] - Hu[n] \conv Res[n]
% where
%     Res[n] = Hc[n] \conv Ht[n]  -  \delta[n]
% and "\conv" denotes convolution.
%    While the above formulation describes the impact of Hu upon
% quantization noise, it is worth noting that the A and B terms also
% can be used to describe the impact of Hu on energy transfer during
% analysis.  Specifically, the amount of receiver (higher frequency)
% subband energy which remains in that band after modulated lifting
% is found by convolving the original receiver subband signal by A[n],
% while the amount of receiver band energy which finds its way into
% the donor (lower frequency) subband during modulated lifting is
% found by applying B[n] to the original receiver band signal.  Thus,
% B[n] determines the amount of aliasing energy which is carried into
% the receiver band, while A[n] describes the amount of energy which
% is left behind in the receiver band.  In an ideal orthonormal system,
% the energies of the A[n] and B[n] filters should sum to 1.  If the
% system were close to orthonormal, we would also expect (at high
% bit-rates) that the quantization noise powers in all subbands should
% be roughly equal (and white), in which case the synthesis objective
% is also to minimize the sum of the energies of the A[n] and B[n]
% filters.  We shall find that this minimization objective cannot
% generally achieve energy sums lower than 1.
%    In summary, our objective is to find Hu, such that
% \sum_n A^2[n] + B^2[n] is minimized.  Equivalently, we seek to
% minimize \int | \hat{A}(\omega) |^2 + | \hat{B}(\omega) |^2  d\omega

Nc = -floor((length(Hc)-1)/2);
Pc = floor(length(Hc)/2);
Nt = -floor(length(Ht)/2);
Pt = floor((length(Ht)-1)/2);
if ((Nc+Nt+Pc+Pt) ~= 0)
   error('Ht and Hc filters have orders which are incompatible with end-to-end zero phase');
end
Nu = Nc;
Pu = Pc;

w = (0:2047)*pi/1024;
Hc_freq = exp(-j*w'*(Nc:Pc))*Hc';
Ht_freq = exp(-j*w'*(Nt:Pt))*Ht';
Res_freq = Hc_freq .* Ht_freq - 1;
desired_freq = conj(Ht_freq) + Hc_freq .* conj(Res_freq);
weight_freq = (Ht_freq .* conj(Ht_freq)) + (Res_freq .* conj(Res_freq));
R_cor = zeros(1,(Pu-Nu+1));
D_vec = zeros(Pu-Nu+1,1);
for n=1:(Pu-Nu+1)
   d_pos = n-1+Nu;
   D_vec(n) = real(exp(j*w*d_pos)*desired_freq) / length(w);   
   r_pos = n-1;
   R_cor(n) = real(exp(j*w*r_pos)*weight_freq) / length(w);
end

R_mat = zeros((Pu-Nu+1),(Pu-Nu+1));
for n=1:(Pu-Nu+1)
   for m=n:(Pu-Nu+1)
      R_mat(m,n) = R_cor(m-n+1);
      R_mat(n,m) = R_cor(m-n+1);
   end
end
Hu = (inv(R_mat)*D_vec)';

Hu_freq = exp(-j*w'*(Nu:Pu))*Hu';

A_freq = 1 - Ht_freq .* Hu_freq;
B_freq = Hc_freq - Res_freq .* Hu_freq;
A_energy = sum(abs(A_freq).^2) / length(w);
B_energy = sum(abs(B_freq).^2) / length(w);
J = A_energy + B_energy;  % Value of objective after optimization
N1_gain = A_energy + sum(abs(Ht_freq).^2)/length(w);
N2_gain = B_energy + sum(abs(Res_freq).^2)/length(w);

Ac_freq = 1 - Ht_freq .* Hc_freq;
Bc_freq = Hc_freq - Res_freq .* Hc_freq;
Ac_energy = sum(abs(Ac_freq).^2) / length(w);
Bc_energy = sum(abs(Bc_freq).^2) / length(w);
Jc = Ac_energy + Bc_energy; % Value of objective, if we use Hc for Hu
N1c_gain = Ac_energy + sum(abs(Ht_freq).^2)/length(w);
N2c_gain = Bc_energy + sum(abs(Res_freq).^2)/length(w);

An_freq = Ac_freq .* 0 + 1;
Bn_freq = Hc_freq;
An_energy = sum(abs(An_freq).^2) / length(w);
Bn_energy = sum(abs(Bn_freq).^2) / length(w);
Jn = An_energy + Bn_energy; % Value of objective, if we use no update at all
N1n_gain = An_energy + sum(abs(Ht_freq).^2)/length(w);
N2n_gain = Bn_energy + sum(abs(Res_freq).^2)/length(w);

% Write comparison statistics
AB_energies = [An_energy  Ac_energy  A_energy;
               Bn_energy  Bc_energy  B_energy;
               Jn  Jc  J]
N1N2_gains = [N1n_gain  N1c_gain  N1_gain;
              N2n_gain  N2c_gain  N2_gain;
              N1n_gain+N2n_gain  N1c_gain+N2c_gain  N1_gain+N2_gain]
           


