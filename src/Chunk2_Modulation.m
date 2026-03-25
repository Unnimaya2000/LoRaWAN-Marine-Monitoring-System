%% =========================================================
%  CHUNK 2 of 5 — CSS Chirp Modulation + Hardware Impairments
%  Generates real baseband IQ chirp signal, then injects:
%    - Carrier Frequency Offset (CFO)
%    - Sampling clock offset
%    - Phase noise
%  Plots: waveform, spectrogram, IQ constellation, PSD
%  Requires: lora_complete.mat  (from Chunk 1)
%  Saves   : lora_complete.mat  (appends chirp signals)
%% =========================================================
clc;
load('lora_complete.mat');
fprintf('--- CHUNK 2: CSS Modulation + Hardware Impairments ---\n\n');

%% -------------------------------------------------------
%  SECTION 2A: Bit Stream → Symbol Mapping
%% -------------------------------------------------------
bits_per_symbol = params.SF;   % SF bits per chirp symbol (10 for SF10)

%% --- STEP 0: LFSR Whitening (Data Scrambling) ---
% XOR the raw bit stream with a pseudo-random sequence from a 9-bit LFSR.
% Polynomial: x^9 + x^5 + 1  (matches LoRa whitening approach)
% Purpose: break up long runs of 0s/1s that degrade FEC and cause spectral lines.
% The same sequence is applied at receiver to de-whiten after decoding.
raw_bits   = bit_stream;
n_raw      = length(raw_bits);

% Generate LFSR sequence of length n_raw
lfsr_state = [1 0 0 1 0 1 0 1 0];  % 9-bit seed (non-zero)
lfsr_seq   = zeros(1, n_raw);
for i = 1:n_raw
    lfsr_seq(i)  = lfsr_state(1);            % output bit
    feedback     = xor(lfsr_state(9), lfsr_state(5));  % taps: 9,5
    lfsr_state   = [feedback, lfsr_state(1:8)];        % shift register
end

% XOR data with whitening sequence
whitened_bits = xor(raw_bits, lfsr_seq);

fprintf('[LFSR Whitener — 9-bit, poly x^9+x^5+1]\n');
fprintf('  Ones before whitening : %d/%d  (%.1f%%)\n', ...
    sum(raw_bits), n_raw, 100*mean(raw_bits));
fprintf('  Ones after  whitening : %d/%d  (%.1f%%)\n', ...
    sum(whitened_bits), n_raw, 100*mean(whitened_bits));

%% --- STEP 1: Hamming(4,8) FEC Encoding ---
% Systematic code: [d1 d2 d3 d4 | p1 p2 p3 p4]
% p1=d1^d2^d3, p2=d1^d2^d4, p3=d1^d3^d4, p4=d2^d3^d4
% Doubles bit count: N bits → 2N bits (code rate = 1/2)
raw_bits  = whitened_bits;   % encode the whitened bits
n_raw     = length(raw_bits);
pad4      = mod(-n_raw, 4);
raw_padded= [raw_bits, zeros(1, pad4)];
n_blocks  = length(raw_padded) / 4;

encoded_bits = zeros(1, n_blocks * 8);
for b = 1:n_blocks
    d  = raw_padded((b-1)*4 + (1:4));
    p1 = xor(xor(d(1),d(2)),d(3));
    p2 = xor(xor(d(1),d(2)),d(4));
    p3 = xor(xor(d(1),d(3)),d(4));
    p4 = xor(xor(d(2),d(3)),d(4));
    encoded_bits((b-1)*8+(1:8)) = [d(1) d(2) d(3) d(4) p1 p2 p3 p4];
end
n_encoded = length(encoded_bits);

fprintf('[Hamming FEC Encoder — CR=4/8]\n');
fprintf('  Input bits    : %d\n', n_raw);
fprintf('  Encoded bits  : %d  (code rate=0.5)\n', n_encoded);

%% --- STEP 2: LoRa Diagonal Bit Interleaver ---
% LoRa spec interleaves bits diagonally across a SF × (4+CR) matrix.
% This ensures that bits from the same FEC codeword are spread across
% different chirp symbols — a burst error corrupting one symbol
% only damages one bit per FEC block, which Hamming can correct.
%
% Matrix layout: SF rows × n_codewords_per_block columns
% Write: column-by-column (FEC codeword bits down each column)
% Read:  diagonal — shift each row by its row index before reading
CR_int   = params.CR;            % coding rate index (1 for 4/5)
SF_int   = params.SF;            % spreading factor (10)
n_cols   = ceil(n_encoded / SF_int);  % number of interleave columns
pad_inter= n_cols * SF_int - n_encoded;
enc_pad  = [encoded_bits, zeros(1, pad_inter)];

% Write column-by-column into SF × n_cols matrix
int_mat  = reshape(enc_pad, SF_int, n_cols);

% Diagonal shift: row r is cyclically shifted left by r positions
for r = 1:SF_int
    int_mat(r,:) = circshift(int_mat(r,:), -(r-1));
end

% Read row-by-row to get interleaved stream
interleaved_bits = reshape(int_mat, 1, []);

fprintf('[LoRa Diagonal Interleaver]\n');
fprintf('  Matrix        : %d rows (SF) x %d cols\n', SF_int, n_cols);
fprintf('  Each row r shifted left by r-1 positions\n');

%% --- STEP 3: Map interleaved bits → symbol values ---
pad_sym   = mod(-length(interleaved_bits), bits_per_symbol);
sym_bits  = [interleaved_bits, zeros(1, pad_sym)];
n_symbols = length(sym_bits) / bits_per_symbol;

symbols_raw = zeros(1, n_symbols);
for s = 1:n_symbols
    idx_s         = (s-1)*bits_per_symbol + (1:bits_per_symbol);
    symbols_raw(s)= bi2de(sym_bits(idx_s), 'left-msb');
end

%% --- STEP 4: Gray Coding ---
% Gray(k) = k XOR floor(k/2)
% Adjacent symbols differ by exactly 1 bit → one-bin FFT error = 1 bit error
symbols = bitxor(symbols_raw, bitshift(symbols_raw, -1));

fprintf('[Gray Coder]\n');
fprintf('  Total symbols : %d  (raw: %d, after FEC+interleave)\n', ...
    n_symbols, ceil(n_raw/bits_per_symbol));
fprintf('  Sample: raw=%d  gray=%d\n', symbols_raw(1), symbols(1));

n_raw_bits     = n_raw;
n_encoded_bits = n_encoded;

%% -------------------------------------------------------
%  SECTION 2B: Generate CSS Chirp Signal (Real IQ Samples)
%  Core LoRa chirp formula:
%    Symbol k shifts the starting frequency bin.
%    phase[n] = 2π × (freq_bin(n) × n) / N
%    where freq_bin(n) = mod(k + n, N)
%% -------------------------------------------------------
fprintf('\n[Generating CSS Chirp IQ Signal...]\n');

chirp_signal = zeros(1, n_symbols * params.N);

for s = 1:n_symbols
    k      = symbols(s);       % symbol value 0..N-1
    % For each chip n: instantaneous frequency bin = mod(k+n, N)
    % Phase = cumulative integral of instantaneous frequency
    n_vec  = (0:params.N-1);
    freq_bins = mod(k + n_vec, params.N);         % circular shift by k
    % Quadratic phase: each chip's phase accumulates from previous
    % phase(n) = 2π/N * sum_{i=0}^{n} freq_bins(i)
    cum_phase  = 2*pi/params.N * cumsum(freq_bins);
    chirp_sym  = exp(1j * cum_phase);
    chirp_signal((s-1)*params.N + (1:params.N)) = chirp_sym;
end

fprintf('  Symbols generated : %d\n', n_symbols);
fprintf('  Total IQ samples  : %d\n', length(chirp_signal));
fprintf('  Signal duration   : %.3f ms\n', length(chirp_signal)/params.BW*1e3);
fprintf('  Mean power        : %.4f (should be ~1.0)\n', mean(abs(chirp_signal).^2));

%% -------------------------------------------------------
%  SECTION 2C: Apply Hardware Impairments (From Code A)
%% -------------------------------------------------------
fprintf('\n[Applying Hardware Impairments]\n');

% --- 2C-i: Carrier Frequency Offset (CFO) ---
% Multiply signal by complex exponential: e^(j*2π*CFO*t)
t_vec   = (0:length(chirp_signal)-1) / params.BW;  % time axis (s)
cfo_rot = exp(1j * 2*pi * hw.CFO_Hz * t_vec);
chirp_cfo = chirp_signal .* cfo_rot;

cfo_norm = hw.CFO_Hz / (params.BW / params.N);
fprintf('  CFO applied      : ±%.1f Hz  (%.4f of symbol BW)\n', hw.CFO_Hz, cfo_norm);

% --- 2C-ii: Sampling Clock Offset ---
% Fractional sample drift — accumulates over packet
% Model: resample signal at slightly different rate
clk_err   = hw.clk_ppm * 1e-6;
n_samples = length(chirp_cfo);
% Create resampled time axis with clock error
t_ideal   = 0:n_samples-1;
t_clk_err = t_ideal * (1 + clk_err);
% Interpolate — linear for efficiency
t_interp  = min(t_clk_err, n_samples-1);  % clip at boundary
xi        = floor(t_interp) + 1;
xi        = min(xi, n_samples);
frac      = t_interp - floor(t_interp);
xi2       = min(xi+1, n_samples);
chirp_clk = chirp_cfo(xi) .* (1-frac) + chirp_cfo(xi2) .* frac;

accum_err = clk_err * length(chirp_clk);
fprintf('  Clock offset     : %d ppm  (%.4f samples drift over packet)\n', ...
    hw.clk_ppm, accum_err);

% --- 2C-iii: Phase Noise (Wiener process model) ---
% Phase increments are iid Gaussian with variance = 2π*linewidth/Fs
pn_sigma  = sqrt(2*pi * hw.pn_linewidth / params.BW);
pn_incr   = pn_sigma * randn(1, length(chirp_clk));
pn_phase  = cumsum(pn_incr);   % Brownian motion in phase
chirp_hw  = chirp_clk .* exp(1j * pn_phase);

pn_penalty_dB = 10*log10(1 + hw.pn_linewidth / (params.BW/params.N));
fprintf('  Phase noise      : linewidth=%.0fHz  penalty=%.4f dB\n', ...
    hw.pn_linewidth, pn_penalty_dB);

% --- 2C-iv: Rx Filter Insertion Loss ---
filter_gain = 10^(-hw.rx_filter_loss_dB/20);
chirp_hw    = chirp_hw * filter_gain;
fprintf('  Rx filter loss   : %.1f dB applied\n', hw.rx_filter_loss_dB);
fprintf('  Final IQ samples : %d  (power: %.4f)\n', ...
    length(chirp_hw), mean(abs(chirp_hw).^2));

%% -------------------------------------------------------
%  SECTION 2D: Reference Dechirp (needed for demodulator)
%  The demodulator multiplies received signal by conjugate
%  of the base up-chirp (k=0), then takes FFT.
%% -------------------------------------------------------
n_ref     = (0:params.N-1);
freq_ref  = mod(0 + n_ref, params.N);
cum_ref   = 2*pi/params.N * cumsum(freq_ref);
base_chirp= exp(1j * cum_ref);
ref_dechirp = conj(base_chirp);   % conjugate = down-chirp reference

fprintf('\n[Reference Dechirp Built]\n');
fprintf('  Base chirp power : %.4f\n', mean(abs(base_chirp).^2));

%% -------------------------------------------------------
%  SECTION 2E: Verify Modulation — Quick Decode at High SNR
%% -------------------------------------------------------
fprintf('\n[Quick Decode Verification at SNR = +20 dB]\n');
snr_test  = 20;
ns_test   = (1/sqrt(2*10^(snr_test/10))) * ...
            (randn(size(chirp_signal)) + 1j*randn(size(chirp_signal)));
rx_test   = chirp_signal + ns_test;

test_ok   = 0;
n_check   = min(8, n_symbols);
sym_check = zeros(1, n_check);
for s = 1:n_check
    seg    = rx_test((s-1)*params.N + (1:params.N));
    dc     = seg .* ref_dechirp;
    [~,pk] = max(abs(fft(dc)));
    sym_check(s) = mod(pk-1, params.N);
    if sym_check(s) == symbols(s); test_ok = test_ok + 1; end
end
fprintf('  Sent    : '); fprintf('%5d ', symbols(1:n_check));   fprintf('\n');
fprintf('  Decoded : '); fprintf('%5d ', sym_check);             fprintf('\n');
fprintf('  Match   : %d / %d\n', test_ok, n_check);

%% -------------------------------------------------------
%  SECTION 2F: Plots
%% -------------------------------------------------------

% ---- Figure 1: Chirp waveform ----
figure('Name','Chunk 2a — Chirp Waveform','Position',[50 50 1100 450]);

subplot(1,3,1);
n_plot = min(2*params.N, length(chirp_signal));
t_us   = (0:n_plot-1) / params.BW * 1e6;
plot(t_us, real(chirp_signal(1:n_plot)), 'b-',  'LineWidth', 1.2); hold on;
plot(t_us, imag(chirp_signal(1:n_plot)), 'r--', 'LineWidth', 1.2);
xlabel('Time (µs)'); ylabel('Amplitude');
title('CSS Chirp Signal — First 2 Symbols');
legend('I (Real)','Q (Imag)','Location','best'); grid on;

subplot(1,3,2);
% Instantaneous frequency of first symbol
first_sym   = chirp_signal(1:params.N);
inst_phase  = unwrap(angle(first_sym));
inst_freq   = diff(inst_phase) / (2*pi) * params.BW / 1e3;  % kHz
t_chip_us   = (0:params.N-2) / params.BW * 1e6;
plot(t_chip_us, inst_freq, 'Color',[0.1 0.5 0.8], 'LineWidth', 2);
xlabel('Time (µs)'); ylabel('Inst. Frequency (kHz)');
title(sprintf('Instantaneous Frequency — Symbol %d (value=%d)', 1, symbols(1)));
grid on; ylim([0 params.BW/1e3]);

subplot(1,3,3);
% Compare clean vs impaired signal IQ
n_iq = min(300, length(chirp_signal));
scatter(real(chirp_signal(1:n_iq)), imag(chirp_signal(1:n_iq)), ...
    10, 'b', 'filled', 'MarkerFaceAlpha',0.4); hold on;
scatter(real(chirp_hw(1:n_iq)), imag(chirp_hw(1:n_iq)), ...
    10, 'r', 'filled', 'MarkerFaceAlpha',0.3);
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
title('IQ Constellation — Clean vs Hardware-Impaired');
legend('Clean TX','After CFO+Phase noise','Location','best');
axis equal; grid on;

sgtitle('LoRa CSS Chirp Signal | SF10 | 125 kHz','FontWeight','bold');

% ---- Figure 2: Spectrogram + PSD ----
figure('Name','Chunk 2b — Spectrogram & PSD','Position',[50 550 1100 400]);

subplot(1,2,1);
n_spec  = min(8*params.N, length(chirp_signal));
win_len = min(params.N, 128);
[~,f_s,t_s,P_s] = spectrogram(chirp_signal(1:n_spec), win_len, ...
    floor(win_len*0.75), 512, params.BW);
imagesc(t_s*1e3, f_s/1e3, 10*log10(P_s+eps));
axis xy; colormap('jet');
xlabel('Time (ms)'); ylabel('Frequency (kHz)');
title('Spectrogram — LoRa Chirp Sweep');
colorbar; caxis([-30 0]);
ylim([0, params.BW/2e3]);

subplot(1,2,2);
seg_len  = min(16*params.N, length(chirp_signal));
[pxx, f_psd] = pwelch(chirp_signal(1:seg_len), [], [], [], params.BW);
f_shift  = (f_psd - max(f_psd)/2) / 1e3;
plot(f_shift, 10*log10(fftshift(pxx)), 'b-', 'LineWidth', 1.8);
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('Power Spectral Density — LoRa Signal');
xline(params.BW/2/1e3, 'r--', '+BW/2'); xline(-params.BW/2/1e3, 'r--', '-BW/2');
grid on;

sgtitle('LoRa Spectral Analysis | SF10 | 125 kHz BW','FontWeight','bold');

%% -------------------------------------------------------
%  Save
%% -------------------------------------------------------
save('lora_complete.mat', 'params', 'hw', 'channel', 'sensor', ...
     'packet', 'payload_bytes', 'bit_stream', 'crc_val', ...
     'airtime_ms', 'n_payload', 'PL', 'N', 'SF', 'BW', 'CR', ...
     'chirp_signal', 'chirp_hw', 'symbols', 'symbols_raw', 'n_symbols', ...
     'bits_per_symbol', 'ref_dechirp', 'base_chirp', ...
     'encoded_bits', 'interleaved_bits', 'n_raw_bits', 'n_encoded_bits', ...
     'n_cols', 'pad_inter', 'pad_sym', 'pad4', ...
     'lfsr_seq', 'whitened_bits');

fprintf('\n✓ Chunk 2 complete. Saved to lora_complete.mat\n');
fprintf('  → Run Chunk3_Channel.m next\n\n');
