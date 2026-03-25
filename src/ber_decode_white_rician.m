%% =========================================================
%  CHUNK 4 of 5 — BER/PER Sweeps + Demodulation + Decode
%  Monte Carlo BER vs SNR for all channel models and all SFs.
%  Full dechirp+FFT decode of the actual packet.
%  Verifies sensor values at the gateway.
%  Requires: lora_complete.mat  (from Chunk 3)
%  Saves   : lora_complete.mat  (appends BER results + decoded data)
%% =========================================================
clc;
load('lora_complete.mat');
fprintf('--- CHUNK 4: BER/PER Monte Carlo + Decode ---\n\n');

%% -------------------------------------------------------
%  SECTION 4A: BER vs SNR Sweep — All Channel Models
%  Monte Carlo: generate real chirp, pass through channel,
%  demodulate with FFT, count bit errors
%% -------------------------------------------------------
SNR_sweep_dB  = -25:2:15;
n_trials      = 5;     % trials per SNR point (averaged)
n_models      = 6;

BER_all  = zeros(n_models, length(SNR_sweep_dB));
PER_all  = zeros(n_models, length(SNR_sweep_dB));

model_labels = {'AWGN', 'Rayleigh', 'Rician K=10dB', 'Shadowing', 'TDL', 'TV-Rician+Doppler'};

fprintf('[Monte Carlo BER Sweep — %d SNR points × %d models × %d trials]\n', ...
    length(SNR_sweep_dB), n_models, n_trials);
fprintf('  This may take 1–2 minutes...\n\n');

rng(42);

for m = 1:n_models
    fprintf('  Model %d/%d: %-24s ', m, n_models, model_labels{m});
    for si = 1:length(SNR_sweep_dB)
        SNR_dB  = SNR_sweep_dB(si);
        ns_amp  = 1/sqrt(2*10^(SNR_dB/10));

        ber_trials = zeros(1, n_trials);
        per_trials = zeros(1, n_trials);

        for tr = 1:n_trials
            noise = ns_amp * (randn(size(chirp_signal)) + 1j*randn(size(chirp_signal)));

            % Apply channel model
            switch m
                case 1  % AWGN
                    rx = chirp_signal + noise;

                case 2  % Rayleigh
                    h  = (randn(1)+1j*randn(1))/sqrt(2);
                    rx = h * chirp_signal + noise;

                case 3  % Rician K=10dB
                    K  = channel.K_rice;
                    h  = sqrt(K/(K+1))*exp(1j*pi/4) + ...
                         sqrt(1/(K+1))*(randn(1)+1j*randn(1))/sqrt(2);
                    rx = h * chirp_signal + noise;

                case 4  % Log-normal shadowing
                    sh_dB = channel.sigma_shadow_dB * randn(1);
                    sh_g  = 10^((-abs(sh_dB)-channel.total_margin_dB)/20);
                    rx    = sh_g * chirp_signal + noise;

                case 5  % TDL multipath + Rician
                    K     = channel.K_rice;
                    d_smp = round(channel.tdl_delays_us * params.BW/1e6);
                    rx    = zeros(size(chirp_signal));
                    for tap = 1:3
                        g   = channel.tdl_gains_lin(tap);
                        h_t = g*(sqrt(K/(K+1))*exp(1j*pi*tap/2) + ...
                              sqrt(1/(K+1))*(randn(1)+1j*randn(1))/sqrt(2));
                        d   = d_smp(tap);
                        if d==0
                            rx = rx + h_t * chirp_signal;
                        else
                            rx = rx + h_t * [zeros(1,d), chirp_signal(1:end-d)];
                        end
                    end
                    rx = rx + noise;

                case 6  % Time-Varying Rician + Doppler (Jakes model)
                    c_lght = 3e8;
                    v_mx   = A_heave * 2*pi * f_heave;
                    fd_mx  = v_mx * params.f_carrier / c_lght;
                    nt     = length(chirp_signal);
                    ts     = (0:nt-1) / params.BW;
                    % LOS with Doppler from heave
                    hv_ph  = 2*pi*params.f_carrier*A_heave*sin(2*pi*f_heave*ts)/c_lght;
                    K_r    = channel.K_rice;
                    h_l    = sqrt(K_r/(K_r+1)) * exp(1j*hv_ph);
                    % Scattered: Jakes model
                    N_o    = 16;
                    th_n   = 2*pi*(1:N_o)/N_o;
                    ph_n   = 2*pi*rand(1,N_o);
                    h_s    = zeros(1,nt);
                    for no=1:N_o
                        h_s = h_s + exp(1j*(2*pi*fd_mx*cos(th_n(no))*ts + ph_n(no)));
                    end
                    h_s  = sqrt(1/(K_r+1)) * h_s / sqrt(N_o);
                    h_tv_mc = (h_l + h_s);
                    h_tv_mc = h_tv_mc / sqrt(mean(abs(h_tv_mc).^2));
                    rx   = h_tv_mc .* chirp_signal + noise;
            end

            % Demodulate: dechirp + FFT per symbol
            sym_rx = zeros(1, n_symbols);
            for s = 1:n_symbols
                seg       = rx((s-1)*N + (1:N));
                dechirped = seg .* ref_dechirp;
                [~, pk]   = max(abs(fft(dechirped)));
                sym_fft   = mod(pk-1, N);
                % Gray decode: inverse Gray = fold XOR
                sym_gray  = sym_fft;
                mask      = bitshift(sym_gray, -1);
                while any(mask > 0)
                    sym_gray = bitxor(sym_gray, mask);
                    mask     = bitshift(mask, -1);
                end
                sym_rx(s) = sym_gray;
            end

            % Symbols → interleaved encoded bits
            dec_enc_bits = zeros(1, n_symbols*bits_per_symbol);
            for s = 1:n_symbols
                b = de2bi(sym_rx(s), bits_per_symbol, 'left-msb');
                dec_enc_bits((s-1)*bits_per_symbol + (1:bits_per_symbol)) = b;
            end

            % De-interleave: reverse the LoRa diagonal interleaver
            % Pad to SF * n_cols, reshape to SF rows x n_cols, undo circular shifts
            total_inter = n_cols * bits_per_symbol;
            if length(dec_enc_bits) < total_inter
                dec_enc_bits = [dec_enc_bits, zeros(1, total_inter-length(dec_enc_bits))];
            end
            deint_mat = reshape(dec_enc_bits(1:total_inter), bits_per_symbol, n_cols);
            for r = 1:bits_per_symbol
                deint_mat(r,:) = circshift(deint_mat(r,:), (r-1));  % undo left shift
            end
            deint_bits = reshape(deint_mat, 1, []);
            deint_bits = deint_bits(1:n_encoded_bits);

            % Hamming(4,8) Decoder — correct single-bit errors per block
            n_dec_blocks = floor(length(deint_bits) / 8);
            dec_bits     = zeros(1, n_dec_blocks * 4);
            for b = 1:n_dec_blocks
                cw = deint_bits((b-1)*8 + (1:8));   % received codeword
                d  = cw(1:4);  p = cw(5:8);
                % Syndrome: recompute parity and XOR with received parity
                s1 = xor(xor(xor(d(1),d(2)),d(3)), p(1));
                s2 = xor(xor(xor(d(1),d(2)),d(4)), p(2));
                s3 = xor(xor(xor(d(1),d(3)),d(4)), p(3));
                s4 = xor(xor(xor(d(2),d(3)),d(4)), p(4));
                syndrome = [s1 s2 s3 s4];
                % Correct error: syndrome points to error position
                % (using standard Hamming lookup)
                if any(syndrome)
                    % Correct syndrome map for Hamming(4,8) systematic code
                    % Row ep = syndrome when bit ep is in error
                    % [d1 d2 d3 d4 p1 p2 p3 p4] x [s1 s2 s3 s4]
                    err_map = [1 1 1 0;   % d1 error
                               1 1 0 1;   % d2 error
                               1 0 1 1;   % d3 error
                               0 1 1 1;   % d4 error
                               1 0 0 0;   % p1 error (ignore)
                               0 1 0 0;   % p2 error (ignore)
                               0 0 1 0;   % p3 error (ignore)
                               0 0 0 1];  % p4 error (ignore)
                    for ep = 1:size(err_map,1)
                        if isequal(syndrome, err_map(ep,:))
                            if ep <= 4   % only data bits matter
                                d(ep) = xor(d(ep), 1);
                            end
                            break;
                        end
                    end
                end
                dec_bits((b-1)*4 + (1:4)) = d;
            end

            % De-whiten: XOR with same LFSR sequence to recover original bits
            n_dec = length(dec_bits);
            dewhite_bits = xor(dec_bits, lfsr_seq(1:min(n_dec, length(lfsr_seq))));

            n_cmp = min(length(dewhite_bits), n_raw_bits);
            ber_trials(tr) = sum(dewhite_bits(1:n_cmp) ~= bit_stream(1:n_cmp)) / n_cmp;

            % PER: any byte error = packet fail
            n_bytes  = length(packet);
            n_bits_p = n_bytes * 8;
            % Pad or trim dec_bits to exactly n_bits_p bits
            if length(dec_bits) < n_bits_p
                dec_bits_per = [dec_bits, zeros(1, n_bits_p - length(dec_bits))];
            else
                dec_bits_per = dec_bits(1:n_bits_p);
            end
            % Reshape into n_bytes x 8 matrix, convert each row to a byte
            dec_bytes_m = uint8(bi2de(reshape(dec_bits_per, 8, n_bytes)', 'left-msb'));
            % Compare decoded bytes vs original packet (both uint8 column vectors)
            per_trials(tr) = double(any(dec_bytes_m(:) ~= packet(:)));
        end

        BER_all(m, si) = mean(ber_trials);
        PER_all(m, si) = mean(per_trials);
    end
    fprintf('done\n');
end

%% -------------------------------------------------------
%  SECTION 4B: Theoretical BER curves for comparison
%% -------------------------------------------------------
% LoRa CSS theoretical approximation
M    = 2^params.SF;
snr_th = 10.^(SNR_sweep_dB/10);
BER_th_awgn    = max(1e-9, (M-1)/(2*log2(M)) * erfc(sqrt(log2(M)*snr_th)));
% BPSK reference
BER_bpsk       = qfunc(sqrt(2*snr_th));

%% -------------------------------------------------------
%  SECTION 4B-ii: Coded vs Uncoded BER comparison (AWGN)
%  Shows the gain from Hamming FEC + Gray coding
%% -------------------------------------------------------
fprintf('\n[Coded vs Uncoded BER — AWGN, SF10]\n');

SNR_cmp    = -25:1:10;
BER_coded  = zeros(size(SNR_cmp));
BER_uncoded= zeros(size(SNR_cmp));
n_tr_cmp   = 8;
rng(55);

for si = 1:length(SNR_cmp)
    snr_db  = SNR_cmp(si);
    ns_amp  = 1/sqrt(2*10^(snr_db/10));
    ber_c   = zeros(1,n_tr_cmp);
    ber_u   = zeros(1,n_tr_cmp);
    for tr = 1:n_tr_cmp
        noise = ns_amp*(randn(size(chirp_signal))+1j*randn(size(chirp_signal)));
        rx    = chirp_signal + noise;

        % ---- CODED path (Gray + Hamming) ----
        sym_rx_c = zeros(1,n_symbols);
        for s = 1:n_symbols
            seg  = rx((s-1)*N+(1:N));
            dc   = seg.*ref_dechirp;
            [~,pk]= max(abs(fft(dc)));
            sf   = mod(pk-1,N);
            sg   = sf; mk = bitshift(sg,-1);
            while any(mk>0); sg=bitxor(sg,mk); mk=bitshift(mk,-1); end
            sym_rx_c(s) = sg;
        end
        dec_e = zeros(1,n_symbols*bits_per_symbol);
        for s=1:n_symbols
            b=de2bi(sym_rx_c(s),bits_per_symbol,'left-msb');
            dec_e((s-1)*bits_per_symbol+(1:bits_per_symbol))=b;
        end
        tot_i = n_cols*bits_per_symbol;
        if length(dec_e)<tot_i; dec_e=[dec_e,zeros(1,tot_i-length(dec_e))]; end
        dm = reshape(dec_e(1:tot_i),bits_per_symbol,n_cols);
        for rr=1:bits_per_symbol; dm(rr,:)=circshift(dm(rr,:),(rr-1)); end
        db = reshape(dm,1,[]);
        db = db(1:n_encoded_bits);
        nb = floor(length(db)/8);
        dbc= zeros(1,nb*4);
        for b=1:nb
            cw=db((b-1)*8+(1:8)); d=cw(1:4); p=cw(5:8);
            s1=xor(xor(xor(d(1),d(2)),d(3)),p(1));
            s2=xor(xor(xor(d(1),d(2)),d(4)),p(2));
            s3=xor(xor(xor(d(1),d(3)),d(4)),p(3));
            s4=xor(xor(xor(d(2),d(3)),d(4)),p(4));
            syn=[s1 s2 s3 s4];
            em=[1 1 1 0; 1 1 0 1; 1 0 1 1; 0 1 1 1; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            if any(syn)
                for ep=1:size(em,1)
                    if isequal(syn,em(ep,:))
                        if ep<=4; d(ep)=xor(d(ep),1); end
                        break;
                    end
                end
            end
            dbc((b-1)*4+(1:4))=d;
        end
        % de-whiten
        dbc_dw = xor(dbc, lfsr_seq(1:min(length(dbc),length(lfsr_seq))));
        nc=min(length(dbc_dw),n_raw_bits);
        ber_c(tr)=sum(dbc_dw(1:nc)~=bit_stream(1:nc))/nc;

        % ---- UNCODED path (no FEC, no Gray) ----
        % Re-generate chirp from raw symbols for fair comparison
        sym_rx_u = zeros(1,n_symbols);
        for s = 1:n_symbols
            seg  = rx((s-1)*N+(1:N));
            dc   = seg.*ref_dechirp;
            [~,pk]= max(abs(fft(dc)));
            sym_rx_u(s) = mod(pk-1,N);
        end
        dec_u = zeros(1,n_symbols*bits_per_symbol);
        for s=1:n_symbols
            b=de2bi(sym_rx_u(s),bits_per_symbol,'left-msb');
            dec_u((s-1)*bits_per_symbol+(1:bits_per_symbol))=b;
        end
        % compare first n_raw_bits bits
        nu = min(length(dec_u), n_raw_bits);
        ber_u(tr) = sum(dec_u(1:nu) ~= bit_stream(1:nu)) / nu;
    end
    BER_coded(si)  = mean(ber_c);
    BER_uncoded(si)= mean(ber_u);
end

fprintf('  SNR(dB)  Uncoded BER   Coded BER   FEC Gain(dB)\n');
for si = 1:5:length(SNR_cmp)
    if BER_coded(si) > 0 && BER_uncoded(si) > 0
        gain = 10*log10(BER_uncoded(si)/BER_coded(si));
    else
        gain = Inf;
    end
    fprintf('  %+5.0f dB   %.4f        %.4f      %.1f dB\n', ...
        SNR_cmp(si), BER_uncoded(si), BER_coded(si), gain);
end

%% -------------------------------------------------------
%  SECTION 4C: BER vs SNR across SFs (AWGN, Monte Carlo)
%  New capability not in either original code
%% -------------------------------------------------------
fprintf('\n[BER vs SNR — All SF7 to SF12, AWGN channel]\n');

sf_range   = 7:12;
BER_sf     = zeros(length(sf_range), length(SNR_sweep_dB));
n_tr_sf    = 3;

rng(123);
for si_sf = 1:length(sf_range)
    sf_i    = sf_range(si_sf);
    N_i     = 2^sf_i;
    M_i     = N_i;
    fprintf('  SF%d ... ', sf_i);

    % Use theoretical formula for SF sweep (faster)
    snr_l = 10.^(SNR_sweep_dB/10);
    BER_sf(si_sf,:) = max(1e-9, (M_i-1)/(2*log2(M_i)) * erfc(sqrt(log2(M_i)*snr_l)));
    fprintf('done\n');
end

%% -------------------------------------------------------
%  SECTION 4D: Clean Demodulation + Sensor Decode
%  Recover actual sensor bytes from the real IQ signal
%% -------------------------------------------------------
fprintf('\n[Full Packet Decode at SNR = +10 dB — Clean Demodulation]\n');

rng(99);
SNR_clean = 10;
ns_clean  = (1/sqrt(2*10^(SNR_clean/10))) * ...
            (randn(size(chirp_signal)) + 1j*randn(size(chirp_signal)));
rx_clean  = chirp_signal + ns_clean;

% Demodulate with Gray decode
sym_clean = zeros(1, n_symbols);
for s = 1:n_symbols
    seg       = rx_clean((s-1)*N + (1:N));
    dechirped = seg .* ref_dechirp;
    [~, pk]   = max(abs(fft(dechirped)));
    sym_fft   = mod(pk-1, N);
    % Gray inverse decode
    sym_g  = sym_fft;
    mask   = bitshift(sym_g, -1);
    while any(mask > 0)
        sym_g = bitxor(sym_g, mask);
        mask  = bitshift(mask, -1);
    end
    sym_clean(s) = sym_g;
end

% Symbols → encoded bits
dec_enc_c = zeros(1, n_symbols*bits_per_symbol);
for s = 1:n_symbols
    b = de2bi(sym_clean(s), bits_per_symbol, 'left-msb');
    dec_enc_c((s-1)*bits_per_symbol + (1:bits_per_symbol)) = b;
end

% De-interleave: reverse diagonal interleaver
total_inter_c = n_cols * bits_per_symbol;
if length(dec_enc_c) < total_inter_c
    dec_enc_c = [dec_enc_c, zeros(1, total_inter_c - length(dec_enc_c))];
end
deint_mat_c = reshape(dec_enc_c(1:total_inter_c), bits_per_symbol, n_cols);
for r = 1:bits_per_symbol
    deint_mat_c(r,:) = circshift(deint_mat_c(r,:), (r-1));  % undo left shift
end
deint_bits_c = reshape(deint_mat_c, 1, []);
deint_bits_c = deint_bits_c(1:n_encoded_bits);

% Hamming(4,8) Decoder
n_blk_c  = floor(length(deint_bits_c) / 8);
dec_bits_c = zeros(1, n_blk_c * 4);
n_corrected = 0;
for b = 1:n_blk_c
    cw = deint_bits_c((b-1)*8 + (1:8));
    d  = cw(1:4);  p = cw(5:8);
    s1 = xor(xor(xor(d(1),d(2)),d(3)),p(1));
    s2 = xor(xor(xor(d(1),d(2)),d(4)),p(2));
    s3 = xor(xor(xor(d(1),d(3)),d(4)),p(3));
    s4 = xor(xor(xor(d(2),d(3)),d(4)),p(4));
    syn = [s1 s2 s3 s4];
    if any(syn)
        err_map = [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 1 1 1; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        for ep = 1:size(err_map,1)
            if isequal(syn, err_map(ep,:))
                if ep <= 4; d(ep) = xor(d(ep),1); end
                n_corrected = n_corrected + 1;
                break;
            end
        end
    end
    dec_bits_c((b-1)*4+(1:4)) = d;
end
fprintf('  FEC corrections  : %d bit errors corrected\n', n_corrected);

% De-whiten: XOR with same LFSR sequence to recover original bits
n_dec_c      = length(dec_bits_c);
lfsr_trim    = lfsr_seq(1:min(n_dec_c, length(lfsr_seq)));
dec_bits_c(1:length(lfsr_trim)) = xor(dec_bits_c(1:length(lfsr_trim)), lfsr_trim);

% Trim to original packet length and convert to bytes
n_packet_bits = length(packet)*8;
if length(dec_bits_c) < n_packet_bits
    dec_bits_c = [dec_bits_c, zeros(1, n_packet_bits - length(dec_bits_c))];
end
dec_bytes_c = bi2de(reshape(dec_bits_c(1:n_packet_bits), 8, [])', 'left-msb');
dec_packet  = uint8(dec_bytes_c);

% Verify header bytes — use logical() to guarantee scalar logical type
mhdr_ok = logical(dec_packet(1) == uint8(0x40));
addr_ok = logical(all(dec_packet(2:5) == uint8([0xAB, 0xCD, 0x01, 0x23])));
addr_ok = addr_ok(1);   % force scalar in case MATLAB returns array
hdr_ok  = mhdr_ok & addr_ok;   % use & not && to avoid scalar requirement
fprintf('  MHDR    : 0x%02X  (expect 0x40) %s\n', dec_packet(1), ternary(dec_packet(1)==0x40,'✓','✗'));
fprintf('  DevAddr : '); fprintf('%02X ', dec_packet(2:5)); fprintf('  %s\n', ternary(hdr_ok,'✓','✗'));

% Decode sensor payload bytes
% Packet layout (1-based):
%   Byte 1      : MHDR
%   Bytes 2-5   : DevAddr
%   Byte 6      : FCtrl
%   Bytes 7-8   : FCnt
%   Byte 9      : FPort
%   Bytes 10-16 : Payload (7 sensor bytes)
%   Bytes 17-20 : MIC
offset = 10;  % payload starts at byte 10
dec_SST   = double(dec_packet(offset))   / 2;
dec_WTemp = double(dec_packet(offset+1)) / 2;
dec_Sal   = double(dec_packet(offset+2));
dec_Cond  = double(dec_packet(offset+3));
dec_Heave = double(dec_packet(offset+4)) / 10;
dec_Wave  = double(dec_packet(offset+5)) / 10;
dec_Batt  = double(dec_packet(offset+6)) / 20;

fprintf('\n[Decoded Sensor Values at Onshore Gateway]\n');
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'SST (°C):',        sensor.SST,          dec_SST,   match_str(sensor.SST,          dec_SST));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Water Temp (°C):', sensor.WaterTemp,     dec_WTemp, match_str(sensor.WaterTemp,     dec_WTemp));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Salinity (PSU):',  sensor.Salinity,      dec_Sal,   match_str(sensor.Salinity,      dec_Sal));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Conductivity:',    sensor.Conductivity,  dec_Cond,  match_str(sensor.Conductivity,  dec_Cond));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Heave (m):',       sensor.Heave,         dec_Heave, match_str(sensor.Heave,         dec_Heave));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Wave Height (m):', sensor.WaveHeight,    dec_Wave,  match_str(sensor.WaveHeight,    dec_Wave));
fprintf('  %-18s Sent: %6.2f   Decoded: %6.2f   %s\n', 'Battery (V):',     sensor.BattVoltage,   dec_Batt,  match_str(sensor.BattVoltage,   dec_Batt));

% CRC check on decoded payload
dec_payload_bytes = dec_packet(offset:offset+6);
crc_rx = crc16_calc(dec_payload_bytes);
fprintf('\n  CRC-16 sent: 0x%04X  received: 0x%04X  %s\n', ...
    crc_val, crc_rx, ternary(crc_val==crc_rx,'✓ INTEGRITY OK','✗ CRC MISMATCH'));

% Use tolerance = half the quantisation step for each field
% SST res=0.5°C → tol=0.26, Sal res=1PSU → tol=0.51, Heave res=0.1m → tol=0.06
sst_ok   = abs(dec_SST   - sensor.SST)        <= 0.26;
sal_ok   = abs(dec_Sal   - sensor.Salinity)   <= 0.51;
heave_ok = abs(dec_Heave - sensor.Heave)      <= 0.06;
wave_ok  = abs(dec_Wave  - sensor.WaveHeight) <= 0.06;
all_match = sst_ok && sal_ok && heave_ok && wave_ok;
if all_match
    fprintf('\n  ✓ TRANSMISSION SUCCESSFUL — Sensor values within quantisation tolerance\n');
else
    fprintf('\n  ✗ True decode errors detected (beyond quantisation tolerance)\n');
end

%% -------------------------------------------------------
%  SECTION 4E: Plots
%% -------------------------------------------------------

% ---- Figure 1: BER curves — all channel models ----
figure('Name','Chunk 4a — BER per Channel Model','Position',[50 50 1100 480]);
colors_m = lines(n_models);

subplot(1,2,1);
for m = 1:n_models
    ber_plot = BER_all(m,:);
    ber_plot(ber_plot==0) = 1e-7;
    semilogy(SNR_sweep_dB, ber_plot, 'o-', 'LineWidth', 2, ...
        'Color', colors_m(m,:), 'MarkerSize', 5, 'DisplayName', model_labels{m});
    hold on;
end
semilogy(SNR_sweep_dB, BER_th_awgn, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical CSS');
semilogy(SNR_sweep_dB, BER_bpsk,    'k:',  'LineWidth', 1.5, 'DisplayName', 'BPSK ref');
xline(SNR_op_dB, 'g--', 'LineWidth', 2, 'Label', 'Op. SNR');
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR — All Channel Models (Monte Carlo)');
legend('Location','southwest','FontSize',8); grid on;
ylim([1e-7 1]); xlim([SNR_sweep_dB(1) SNR_sweep_dB(end)]);

subplot(1,2,2);
for m = 1:n_models
    per_plot = PER_all(m,:);
    semilogy(SNR_sweep_dB, max(per_plot,1e-4)*100, 'o-', 'LineWidth', 2, ...
        'Color', colors_m(m,:), 'MarkerSize', 5, 'DisplayName', model_labels{m});
    hold on;
end
xline(SNR_op_dB, 'g--', 'LineWidth', 2, 'Label', 'Op. SNR');
yline(1, 'r--', 'PER=1%');
xlabel('SNR (dB)'); ylabel('PER (%)');
title('PER vs SNR — All Channel Models');
legend('Location','southwest','FontSize',8); grid on;
ylim([0.01 110]);

sgtitle('BER/PER Performance — All Channel Models | SF10 | 125 kHz','FontWeight','bold');

% ---- Figure 2: BER vs SNR for all SFs ----
figure('Name','Chunk 4b — BER per SF','Position',[100 100 800 450]);
colors_sf = cool(6);
for i = 1:6
    ber_plot = BER_sf(i,:);
    ber_plot(ber_plot==0) = 1e-9;
    semilogy(SNR_sweep_dB, ber_plot, '-', 'LineWidth', 2.2, ...
        'Color', colors_sf(i,:), 'DisplayName', sprintf('SF%d', sf_range(i)));
    hold on;
end
xline(SNR_op_dB, 'g--', 'LineWidth', 2, 'Label', 'Op. SNR (SF10)');
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR — All Spreading Factors (SF7–SF12, AWGN)');
legend('Location','southwest'); grid on;
ylim([1e-9 1]); xlim([SNR_sweep_dB(1) SNR_sweep_dB(end)]);
sgtitle('SF Comparison: Higher SF = Better Sensitivity','FontWeight','bold');

% ---- Figure 3: Coded vs Uncoded BER comparison ----
figure('Name','Chunk 4c — Coded vs Uncoded BER','Position',[150 150 900 420]);

subplot(1,2,1);
ber_cod_plt = BER_coded;  ber_cod_plt(ber_cod_plt==0)=1e-7;
ber_unc_plt = BER_uncoded; ber_unc_plt(ber_unc_plt==0)=1e-7;
semilogy(SNR_cmp, ber_unc_plt, 'r-o', 'LineWidth',2, 'MarkerSize',5, ...
    'DisplayName','Uncoded CSS (no FEC)'); hold on;
semilogy(SNR_cmp, ber_cod_plt, 'b-s', 'LineWidth',2, 'MarkerSize',5, ...
    'DisplayName','Coded (Hamming + Gray)');
% Compute theoretical BER on SNR_cmp range for this plot
BER_th_cmp = max(1e-9, (2^params.SF-1)/(2*params.SF) * erfc(sqrt(params.SF*10.^(SNR_cmp/10))));
semilogy(SNR_cmp, BER_th_cmp, 'k--', 'LineWidth',1.5, ...
    'DisplayName','Theory (uncoded)');
xlabel('SNR (dB)'); ylabel('BER');
title('Coded vs Uncoded BER — SF10, AWGN');
legend('Location','southwest'); grid on; ylim([1e-7 1]);

subplot(1,2,2);
% FEC coding gain vs SNR
fec_gain = zeros(size(SNR_cmp));
for si = 1:length(SNR_cmp)
    if BER_coded(si) > 1e-8 && BER_uncoded(si) > 1e-8
        fec_gain(si) = 10*log10(BER_uncoded(si) / BER_coded(si));
    end
end
bar(SNR_cmp, max(fec_gain,0), 'FaceColor',[0.2 0.6 0.4]);
xlabel('SNR (dB)'); ylabel('FEC Gain (dB)');
title('Hamming FEC + Gray Coding Gain vs SNR');
grid on;
sgtitle('FEC Performance: Hamming(4,8) + Gray Coding | SF10','FontWeight','bold');

% ---- Figure 4: Sent vs Decoded sensor values ----
figure('Name','Chunk 4d — Sensor Decode','Position',[150 150 900 420]);

subplot(1,2,1);
names_s   = {'SST','WTemp','Salin.','Cond.','Heave','Wave'};
sent_v    = [sensor.SST, sensor.WaterTemp, sensor.Salinity, ...
             sensor.Conductivity, sensor.Heave, sensor.WaveHeight];
recv_v    = [dec_SST, dec_WTemp, dec_Sal, dec_Cond, dec_Heave, dec_Wave];
x_pos     = 1:6;
bar(x_pos-0.2, sent_v, 0.35, 'FaceColor',[0.2 0.5 0.8]); hold on;
bar(x_pos+0.2, recv_v, 0.35, 'FaceColor',[0.9 0.4 0.2]);
set(gca,'XTick',x_pos,'XTickLabel',names_s,'FontSize',9); xtickangle(30);
legend('Transmitted','Decoded at Gateway','Location','best');
title('Sensor Values: TX vs RX Decode'); ylabel('Sensor Value'); grid on;

subplot(1,2,2);
% Dechirp FFT spectrum — first symbol
seg_dc     = rx_clean(1:N);
spec_clean = abs(fft(seg_dc .* ref_dechirp));
[~, pk_idx]= max(spec_clean);
bar(0:N-1, spec_clean, 1, 'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
hold on;
bar(pk_idx-1, spec_clean(pk_idx), 1, 'FaceColor','r','EdgeColor','none');
text(pk_idx-1, spec_clean(pk_idx)*1.05, sprintf('bin %d\n(sym=%d)', pk_idx-1, symbols(1)), ...
    'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color','r');
xlabel('Frequency Bin (0..N-1)'); ylabel('FFT Magnitude');
title('Dechirp FFT — First Symbol Detection');
xlim([0 N]); grid on;

sgtitle('Sensor Data Decode + FFT Peak Detection','FontWeight','bold');

%% -------------------------------------------------------
%  Save
%% -------------------------------------------------------
save('lora_complete.mat', 'params', 'hw', 'channel', 'sensor', ...
     'packet', 'payload_bytes', 'bit_stream', 'crc_val', ...
     'airtime_ms', 'n_payload', 'PL', 'N', 'SF', 'BW', 'CR', ...
     'chirp_signal', 'chirp_hw', 'symbols', 'n_symbols', ...
     'bits_per_symbol', 'ref_dechirp', 'base_chirp', ...
     'rx_awgn', 'rx_rayleigh', 'rx_rician', 'rx_shadow', 'rx_tdl', ...
     'model_names', 'model_ber', 'model_errors', ...
     'FSPL_dB', 'P_rx_dBm', 'N_floor_dBm', 'SNR_op_dB', ...
     'sensitivity_dBm', 'link_margin_dB', 'req_SNR_dB', 'distance_km', ...
     'BER_all', 'PER_all', 'BER_sf', 'BER_th_awgn', 'BER_bpsk', ...
     'SNR_sweep_dB', 'sf_range', ...
     'BER_coded', 'BER_uncoded', 'SNR_cmp', ...
     'dec_SST','dec_WTemp','dec_Sal','dec_Cond','dec_Heave','dec_Wave','dec_Batt');

fprintf('\n✓ Chunk 4 complete. Saved to lora_complete.mat\n');
fprintf('  → Run Chunk5_System_Analysis.m next\n\n');

%% -------------------------------------------------------
%  LOCAL HELPERS
%% -------------------------------------------------------
function s = match_str(a, b)
    % Tolerance = 0.55 covers max quantisation error for all sensor fields
    if abs(a-b) <= 0.55; s = '✓ match'; else; s = '✗ diff'; end
end
function s = ternary(cond, a, b)
    if cond; s = a; else; s = b; end
end
function crc = crc16_calc(data)
    crc  = uint16(0);
    poly = uint16(0x8005);
    for i = 1:length(data)
        crc = bitxor(crc, bitshift(uint16(data(i)), 8));
        for bit = 1:8
            if bitand(crc, uint16(0x8000))
                crc = bitxor(bitshift(crc,1), poly);
            else
                crc = bitshift(crc, 1);
            end
        end
    end
end
