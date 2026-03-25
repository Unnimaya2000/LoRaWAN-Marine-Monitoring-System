%% =========================================================
%  CHUNK 3 of 5 — Realistic Channel Models
%  Applies all 5 channel models to the real IQ chirp signal:
%    1. AWGN only
%    2. Rayleigh flat fading + AWGN
%    3. Rician fading (K=10dB, sea LOS) + AWGN
%    4. Log-normal shadowing + FSPL margins + AWGN
%    5. TDL 3-tap multipath + Rician + AWGN
%  Requires: lora_complete.mat  (from Chunk 2)
%  Saves   : lora_complete.mat  (appends channel outputs)
%% =========================================================
clc;
load('lora_complete.mat');
fprintf('--- CHUNK 3: Realistic Channel Models ---\n\n');

%% -------------------------------------------------------
%  SECTION 3A: Link Budget & Operating SNR
%% -------------------------------------------------------
distance_km   = 15;
freq_Hz       = params.f_carrier;

% Free Space Path Loss (Friis formula)
FSPL_dB   = 20*log10(distance_km*1e3) + 20*log10(freq_Hz) - 147.55;

% Received power
P_rx_dBm  = params.P_tx_dBm + params.Gtx_dBi + params.Grx_dBi - FSPL_dB;

% Physics-based noise floor (From Code A — more accurate)
N_floor_dBm = 10*log10(params.k_B * params.T_K * params.BW) + 30 + params.NF_dB;

% Required SNR for SF10 LoRa
req_SNR_dB  = -7.5 - 2.5*(params.SF - 7);   % = -17.5 dB for SF10
sensitivity_dBm = N_floor_dBm + req_SNR_dB;

% Actual operating SNR at 15 km
SNR_op_dB   = P_rx_dBm - N_floor_dBm;

% Link margin
link_margin_dB = P_rx_dBm - sensitivity_dBm - channel.total_margin_dB;

fprintf('[Link Budget — 15 km Arabian Sea Path]\n');
fprintf('  Tx Power         : %d dBm\n',     params.P_tx_dBm);
fprintf('  Tx Antenna       : +%.2f dBi\n',  params.Gtx_dBi);
fprintf('  Rx Antenna       : +%.2f dBi\n',  params.Grx_dBi);
fprintf('  FSPL @ 15 km     : %.2f dB\n',    FSPL_dB);
fprintf('  Rx Power         : %.2f dBm\n',   P_rx_dBm);
fprintf('  Noise Floor      : %.2f dBm\n',   N_floor_dBm);
fprintf('  Operating SNR    : %.2f dB\n',    SNR_op_dB);
fprintf('  Sensitivity      : %.2f dBm  (req SNR = %.1f dB)\n', sensitivity_dBm, req_SNR_dB);
fprintf('  Practical margins: %.1f dB  (rain+sea+tilt)\n', channel.total_margin_dB);
fprintf('  Link Margin      : %.2f dB  ', link_margin_dB);
if link_margin_dB > 0
    fprintf('✓ LINK OK\n\n');
else
    fprintf('✗ LINK MARGINAL\n\n');
end

%% -------------------------------------------------------
%  SECTION 3B: Channel Model Functions
%  All models take the clean chirp + SNR_dB and return
%  the received signal ready for demodulation
%% -------------------------------------------------------

% Working SNR for channel comparison (use operating SNR)
SNR_dB_ch = SNR_op_dB;
rng(42);

fprintf('[Applying Channel Models at SNR = %.1f dB]\n\n', SNR_dB_ch);

% Noise generator utility
noise_gen = @(sig, snr_db) ...
    (1/sqrt(2*10^(snr_db/10))) * ...
    (randn(size(sig)) + 1j*randn(size(sig)));

% NOTE: Channel models applied to chirp_signal (clean IQ).
% chirp_hw (with CFO/phase noise) is saved separately for waveform analysis.
% Using clean signal ensures demodulator can recover symbols correctly.

% ---- MODEL 1: AWGN only ----
noise_awgn    = noise_gen(chirp_signal, SNR_dB_ch);
rx_awgn       = chirp_signal + noise_awgn;
fprintf('  [1] AWGN only applied\n');

% ---- MODEL 2: Rayleigh flat fading + AWGN ----
h_ray         = (randn(1) + 1j*randn(1)) / sqrt(2);  % CN(0,1)
noise_ray     = noise_gen(chirp_signal, SNR_dB_ch);
rx_rayleigh   = h_ray * chirp_signal + noise_ray;
fprintf('  [2] Rayleigh fading: |h|=%.4f, angle=%.2f rad\n', ...
    abs(h_ray), angle(h_ray));

% ---- MODEL 3: Rician fading + AWGN ----
K_rice        = channel.K_rice;
h_los         = sqrt(K_rice/(K_rice+1)) * exp(1j*pi/4);
h_scatter     = sqrt(1/(K_rice+1)) * (randn(1)+1j*randn(1))/sqrt(2);
h_rician      = h_los + h_scatter;
noise_ric     = noise_gen(chirp_signal, SNR_dB_ch);
rx_rician     = h_rician * chirp_signal + noise_ric;
fprintf('  [3] Rician fading:   K=%.0f dB, |h_LOS|=%.4f, |h_total|=%.4f\n', ...
    channel.K_rice_dB, abs(h_los), abs(h_rician));

% ---- MODEL 4: Log-normal shadowing + AWGN ----
shadow_dB     = channel.sigma_shadow_dB * randn(1);
shadow_gain   = 10^((-abs(shadow_dB) - channel.total_margin_dB)/20);
noise_shad    = noise_gen(chirp_signal, SNR_dB_ch);
rx_shadow     = shadow_gain * chirp_signal + noise_shad;
fprintf('  [4] Shadowing:       draw=%.2f dB, gain=%.4f\n', shadow_dB, shadow_gain);

% ---- MODEL 5: TDL 3-tap multipath + Rician + AWGN ----
delay_samples = round(channel.tdl_delays_us * params.BW / 1e6);
rx_tdl        = zeros(size(chirp_signal));

for tap = 1:length(delay_samples)
    d     = delay_samples(tap);
    g_lin = channel.tdl_gains_lin(tap);
    h_tap_los     = sqrt(K_rice/(K_rice+1)) * exp(1j*2*pi*tap/3);
    h_tap_scatter = sqrt(1/(K_rice+1)) * (randn(1)+1j*randn(1))/sqrt(2);
    h_tap         = g_lin * (h_tap_los + h_tap_scatter);
    if d == 0
        rx_tdl = rx_tdl + h_tap * chirp_signal;
    else
        sig_delayed = [zeros(1,d), chirp_signal(1:end-d)];
        rx_tdl      = rx_tdl + h_tap * sig_delayed;
    end
end
noise_tdl   = noise_gen(chirp_signal, SNR_dB_ch);
rx_tdl      = rx_tdl + noise_tdl;
fprintf('  [5] TDL multipath:   %d taps, delays=[%s] samples\n', ...
    length(delay_samples), num2str(delay_samples));

% Check ISI risk
max_delay_us = max(channel.tdl_delays_us);
fprintf('      Max delay=%.1f µs vs symbol=%.1f ms → ISI ratio=%.6f (negligible)\n', ...
    max_delay_us, params.T_symbol*1e3, (max_delay_us*1e-6)/params.T_symbol);

% ---- MODEL 6: Time-Varying Rician + Doppler (Heaving Buoy) ----
% The buoy heaves vertically with sinusoidal motion.
% This creates a time-varying Doppler shift on the LOS path.
%
% Physics:
%   Heave motion: z(t) = A_heave * sin(2*pi*f_heave*t)
%   Heave velocity: v(t) = A_heave * 2*pi*f_heave * cos(2*pi*f_heave*t)
%   Instantaneous Doppler: f_d(t) = v(t) * f_carrier / c
%   Maximum Doppler: f_d_max = A_heave * 2*pi*f_heave * f_carrier / c
%
% Jake's model for scattered component:
%   h_scatter(t) = (1/sqrt(N_osc)) * sum_n cos(2*pi*f_d_max*cos(theta_n)*t + phi_n)
%   where theta_n = 2*pi*n/N_osc, phi_n = random phase
%   This generates a correlated complex Gaussian process with the
%   correct Doppler power spectrum (U-shaped/Clarke's spectrum)

c_light     = 3e8;
f_heave     = 0.1;          % wave frequency (Hz) — 10s wave period
A_heave     = sensor.Heave; % heave amplitude = 1.3 m (from sensor)
v_max       = A_heave * 2*pi * f_heave;   % max heave velocity (m/s)
f_d_max     = v_max * params.f_carrier / c_light;   % max Doppler (Hz)
T_coherence = 0.423 / f_d_max;           % Clarke coherence time (s)

fprintf('\n  [6] Time-Varying Rician + Doppler (heaving buoy)\n');
fprintf('      Heave: A=%.2fm at %.2fHz  →  v_max=%.3f m/s\n', ...
    A_heave, f_heave, v_max);
fprintf('      Max Doppler shift : %.4f Hz\n', f_d_max);
fprintf('      Coherence time    : %.1f ms  (symbol=%.1f ms)\n', ...
    T_coherence*1e3, params.T_symbol*1e3);
fprintf('      Channel changes   : %.2f times across packet\n', ...
    n_symbols * params.T_symbol / T_coherence);

% Time axis for the entire packet (sample-by-sample)
n_total = length(chirp_signal);
t_samp  = (0:n_total-1) / params.BW;   % time per sample (s)

% --- LOS component: Doppler from sinusoidal heave ---
% Phase from heave displacement: phi(t) = 2*pi*f_carrier * z(t)/c
%   = 2*pi*f_carrier * A_heave * sin(2*pi*f_heave*t) / c
heave_phase  = 2*pi * params.f_carrier * A_heave * ...
               sin(2*pi * f_heave * t_samp) / c_light;
K_rice       = channel.K_rice;
los_amp      = sqrt(K_rice / (K_rice + 1));
h_los_tv     = los_amp * exp(1j * heave_phase);   % time-varying LOS

% --- Scattered component: Jake's model (N_osc oscillators) ---
N_osc   = 16;   % number of oscillators (higher = smoother spectrum)
rng_state = rng;  % save RNG state for reproducibility
theta_n = 2*pi * (1:N_osc) / N_osc;   % angle of arrival per oscillator
phi_n   = 2*pi * rand(1, N_osc);       % random initial phases

scat_amp  = sqrt(1 / (K_rice + 1));
h_scat_tv = zeros(1, n_total);
for n_osc = 1:N_osc
    f_osc        = f_d_max * cos(theta_n(n_osc));  % Doppler for this oscillator
    h_scat_tv    = h_scat_tv + ...
        exp(1j * (2*pi * f_osc * t_samp + phi_n(n_osc)));
end
h_scat_tv = scat_amp * h_scat_tv / sqrt(N_osc);

% Combined time-varying channel
h_tv    = h_los_tv + h_scat_tv;

% Normalise so mean power = 1
h_tv    = h_tv / sqrt(mean(abs(h_tv).^2));

fprintf('      h_tv power (mean) : %.4f (normalised to 1.0)\n', mean(abs(h_tv).^2));
fprintf('      h_tv power std    : %.4f (variation across packet)\n', std(abs(h_tv).^2));

% Apply time-varying channel + AWGN
noise_tv     = noise_gen(chirp_signal, SNR_dB_ch);
rx_tv_rician = h_tv .* chirp_signal + noise_tv;

%% -------------------------------------------------------
%  SECTION 3C: Demodulate all 5 channel outputs
%  Quick symbol decode to measure error count per model
%% -------------------------------------------------------
fprintf('\n[Quick Demodulation Across All Channel Models]\n');

rx_signals   = {rx_awgn, rx_rayleigh, rx_rician, rx_shadow, rx_tdl, rx_tv_rician};
model_names  = {'AWGN', 'Rayleigh', 'Rician K=10dB', 'Shadowing', 'TDL Multipath', 'TV-Rician+Doppler'};
model_errors = zeros(1, 6);
model_ber    = zeros(1, 6);

for m = 1:6
    rx_m   = rx_signals{m};
    n_chk  = min(n_symbols, floor(length(rx_m)/params.N));
    sym_rx = zeros(1, n_chk);
    for s = 1:n_chk
        idx_start = (s-1)*params.N + 1;
        idx_end   = min(s*params.N, length(rx_m));
        seg       = rx_m(idx_start:idx_end);
        if length(seg) < params.N
            seg = [seg, zeros(1, params.N - length(seg))];
        end
        dechirped = seg .* ref_dechirp;
        [~, pk]   = max(abs(fft(dechirped)));
        % Gray decode to get raw symbol value
        sym_fft = mod(pk-1, params.N);
        sym_g   = sym_fft;
        mask    = bitshift(sym_g, -1);
        while any(mask > 0)
            sym_g = bitxor(sym_g, mask);
            mask  = bitshift(mask, -1);
        end
        sym_rx(s) = sym_g;
    end
    % Compare against transmitted symbols_raw (after whitening+interleaving+Gray)
    % Use Gray-decoded symbols vs symbols_raw (pre-Gray, post-interleave)
    n_cmp_sym = min(length(sym_rx), length(symbols_raw));
    sym_errs  = sum(sym_rx(1:n_cmp_sym) ~= symbols_raw(1:n_cmp_sym));
    % SER (symbol error rate) and approximate BER
    ser = sym_errs / n_cmp_sym;
    % BER approximation: for Gray-coded M-ary: BER ≈ SER / log2(M)
    approx_ber = ser / bits_per_symbol;
    model_errors(m) = sym_errs;
    model_ber(m)    = approx_ber;
    fprintf('  %-22s : %d sym errors  SER=%.4f  BER≈%.4f\n', ...
        model_names{m}, sym_errs, ser, approx_ber);
end

%% -------------------------------------------------------
%  SECTION 3D: Plots
%% -------------------------------------------------------

% ---- Figure 1: Received signal power across models ----
figure('Name','Chunk 3 — Channel Models','Position',[50 50 1200 500]);

subplot(2,3,1);
% Channel model power comparison
rx_powers = zeros(1,6);
for m=1:6; rx_powers(m) = 10*log10(mean(abs(rx_signals{m}).^2)+eps); end
bar_colors = [0.3 0.7 0.3; 0.9 0.4 0.2; 0.2 0.5 0.9; 0.8 0.2 0.8; 0.9 0.7 0.1; 0.2 0.8 0.8];
b = bar(1:6, rx_powers, 0.6, 'FaceColor','flat');
for i=1:6; b.CData(i,:)=bar_colors(i,:); end
xlbls = {'AWGN','Rayleigh','Rician','Shadow','TDL','TV-Rician'};
set(gca,'XTick',1:6,'XTickLabel',xlbls,'FontSize',8);
xtickangle(30); ylabel('Mean Rx Power (dB)');
title('Received Signal Power by Channel Model'); grid on;

subplot(2,3,2);
bar(1:6, model_ber*100, 0.6, 'FaceColor','flat', 'CData', bar_colors);
set(gca,'XTick',1:6,'XTickLabel',xlbls,'FontSize',8);
xtickangle(30); ylabel('BER (%)');
title(sprintf('BER per Channel Model  (SNR=%.1f dB)', SNR_dB_ch));
grid on;

% ---- TDL tap power visualisation ----
subplot(2,3,3);
stem(delay_samples, channel.tdl_gains_dB, 'filled', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Delay (samples)'); ylabel('Tap Gain (dB)');
title('TDL Channel — Tap Profile');
for i=1:3
    text(delay_samples(i)+0.1, channel.tdl_gains_dB(i)+0.5, ...
        sprintf('%.0fµs / %.0fdB', channel.tdl_delays_us(i), channel.tdl_gains_dB(i)), ...
        'FontSize',8);
end
grid on;

% ---- Received signal envelopes ----
t_plot = (0:params.N-1)/params.BW*1e6;  % µs for one symbol
for m = 1:2
    subplot(2,3,3+m);
    rx_env = abs(rx_signals{m}(1:params.N));
    plot(t_plot, rx_env, 'LineWidth', 1.2);
    xlabel('Time (µs)'); ylabel('|r(t)|');
    title(sprintf('Signal Envelope — %s', model_names{m})); grid on;
end

subplot(2,3,6);
% Time-varying channel amplitude across entire packet
t_packet_ms = (0:n_total-1) / params.BW * 1e3;
plot(t_packet_ms, abs(h_tv), 'b-', 'LineWidth', 1.2, 'DisplayName','|h(t)| TV-Rician'); hold on;
plot(t_packet_ms, abs(h_los_tv), 'r--', 'LineWidth', 1.0, 'DisplayName','LOS component');
% Mark symbol boundaries
for s_idx = 0:n_symbols
    xline(s_idx*params.T_symbol*1e3, 'k:', 'Alpha', 0.3);
end
xlabel('Time (ms)'); ylabel('|h(t)|');
title(sprintf('Time-Varying Channel Amplitude\n(f_d=%.2fHz, T_c=%.0fms)', f_d_max, T_coherence*1e3));
legend('Location','best'); grid on;

sgtitle('Channel Model Comparison — Arabian Sea','FontWeight','bold');

%% -------------------------------------------------------
%  Save
%% -------------------------------------------------------
save('lora_complete.mat', 'params', 'hw', 'channel', 'sensor', ...
     'packet', 'payload_bytes', 'bit_stream', 'crc_val', ...
     'airtime_ms', 'n_payload', 'PL', 'N', 'SF', 'BW', 'CR', ...
     'chirp_signal', 'chirp_hw', 'symbols', 'n_symbols', ...
     'bits_per_symbol', 'ref_dechirp', 'base_chirp', ...
     'rx_awgn', 'rx_rayleigh', 'rx_rician', 'rx_shadow', 'rx_tdl', 'rx_tv_rician', ...
     'h_tv', 'h_los_tv', 'f_d_max', 'T_coherence', 'f_heave', 'A_heave', ...
     'model_names', 'model_ber', 'model_errors', ...
     'FSPL_dB', 'P_rx_dBm', 'N_floor_dBm', 'SNR_op_dB', ...
     'sensitivity_dBm', 'link_margin_dB', 'req_SNR_dB', 'distance_km');

fprintf('\n✓ Chunk 3 complete. Saved to lora_complete.mat\n');
fprintf('  → Run Chunk4_BER_Decode.m next\n\n');
