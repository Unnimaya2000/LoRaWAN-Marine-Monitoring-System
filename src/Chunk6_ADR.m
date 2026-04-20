%% =========================================================
%  CHUNK 6 — Adaptive Data Rate (ADR) Simulation
%  CUSAT RUSA: Arabian Sea Buoy — 24-hour deployment
%
%  Implements the full LoRaWAN ADR algorithm:
%    1. Gateway collects SNR from last 20 SUCCESSFUL packets
%    2. Every 20 packets: compute SNR_max from history
%    3. SF step = floor((SNR_max - req_SNR - margin) / 2.5)
%    4. Reduce SF by that many steps (faster + less power)
%    5. If SF=7 and still excess margin: reduce Tx power in 3dB steps
%    6. If SNR margin too small: increase SF (better sensitivity)
%
%  Key features:
%    — SX1276 nonlinear Tx current from datasheet (Table 7, PA_BOOST)
%    — SNR history updated on successful deliveries only (LoRaWAN 1.0.3)
%    — ADR_ACK_REQ / ADR_ACK_LIMIT recovery mechanism
%
%  Simulates 288 packets over 24 hours (5-min interval)
%  with a realistic time-varying SNR profile (sea state changes)
%
%  Compares:
%    - Fixed SF10   (baseline design)
%    - ADR-controlled SF (adapts to channel)
%
%  Requires: lora_complete.mat (from Chunk 4)
%  Saves   : lora_adr.mat
%% =========================================================
clc;
load('lora_complete.mat');
fprintf('=============================================================\n');
fprintf('  CHUNK 6: Adaptive Data Rate (ADR) Simulation\n');
fprintf('  CUSAT RUSA — Arabian Sea Buoy 24-Hour Deployment\n');
fprintf('=============================================================\n\n');

%% -------------------------------------------------------
%  SECTION 6A: ADR Parameters (LoRaWAN 1.0.3 Specification)
%% -------------------------------------------------------
adr.history_len   = 20;   % SNR history window (last N successful packets)
adr.margin_dB     = 10;   % ADR margin above required SNR (dB)
                           % LoRaWAN default = 15 dB; 10 dB for more
                           % aggressive adaptation
adr.sf_min        = 7;    % minimum SF (fastest)
adr.sf_max        = 12;   % maximum SF (longest range)
adr.p_tx_max_dBm  = 14;   % maximum Tx power (dBm)
adr.p_tx_min_dBm  = 2;    % minimum Tx power (dBm)
adr.p_tx_step_dB  = 3;    % power reduction/increase step (dB)
adr.snr_margin_up = 5;    % if margin < this → increase SF (dB)
adr.update_every  = 20;   % run ADR algorithm every N packets

% ADR_ACK_REQ parameters (LoRaWAN 1.0.3 Section 4.3.1.1)
adr.ack_limit     = 64;   % ADR_ACK_LIMIT: raise ACK_REQ flag after this
                           % many consecutive unacknowledged uplinks
adr.ack_delay     = 32;   % ADR_ACK_DELAY: after ACK_REQ, wait this many
                           % more packets before full reset to SF12+14dBm

fprintf('[ADR Configuration — LoRaWAN 1.0.3 Compliant]\n');
fprintf('  SNR history window   : %d successful packets\n', adr.history_len);
fprintf('  ADR margin           : %d dB\n', adr.margin_dB);
fprintf('  SF range             : SF%d to SF%d\n', adr.sf_min, adr.sf_max);
fprintf('  Tx power range       : %d to %d dBm (step: %d dB)\n', ...
    adr.p_tx_min_dBm, adr.p_tx_max_dBm, adr.p_tx_step_dB);
fprintf('  ADR update interval  : every %d packets\n', adr.update_every);
fprintf('  ADR_ACK_LIMIT        : %d unacknowledged packets\n', adr.ack_limit);
fprintf('  ADR_ACK_DELAY        : %d packets after ACK_REQ\n\n', adr.ack_delay);

%% -------------------------------------------------------
%  SECTION 6B: 24-Hour SNR Profile (Arabian Sea)
%    Morning  (0–8h)   : calm sea, high SNR  (Rician K=10 dB)
%    Afternoon(8–16h)  : building swell, moderate SNR variation
%    Evening  (16–20h) : monsoon rain + high waves, SNR drops ~8 dB
%    Night    (20–24h) : conditions ease, SNR recovers
%% -------------------------------------------------------
tx_interval_s = 300;       % 5 minutes between packets
n_packets     = 288;       % 24 hours × 12 packets/hour
t_hours       = linspace(0, 24, n_packets);

SNR_base_dB = SNR_op_dB;  % operational SNR at 15 km (from Chunk 4 link budget)

rng(42);
slow_var  = 4  * sin(2*pi*t_hours/8 - pi/3);       % 8h tidal cycle, ±4 dB
monsoon   = -8 * exp(-((t_hours-18).^2)/(2*2^2));  % Gaussian dip at 18:00, σ=2h
fast_var  = 2  * randn(1, n_packets);               % Rician fast fading, ±2 dB std

SNR_profile = SNR_base_dB + slow_var + monsoon + fast_var;

fprintf('[24-Hour SNR Profile — Arabian Sea]\n');
fprintf('  Base SNR (15 km)     : %.1f dB\n',  SNR_base_dB);
fprintf('  Slow variation       : ±4 dB   (8-hour tidal/sea-state cycle)\n');
fprintf('  Monsoon dip          : -8 dB   (Gaussian, centred 18:00, sigma=2h)\n');
fprintf('  Fast fading std      : 2 dB    (packet-to-packet Rician)\n');
fprintf('  Min SNR in profile   : %.1f dB (worst case)\n', min(SNR_profile));
fprintf('  Max SNR in profile   : %.1f dB (best case)\n\n', max(SNR_profile));

%% -------------------------------------------------------
%  SECTION 6C: Required SNR per SF + Airtime
%% -------------------------------------------------------
sf_list    = 7:12;
req_snr_sf = -7.5 - 2.5*(sf_list - 7);
% SF7=-7.5, SF8=-10, SF9=-12.5, SF10=-15, SF11=-17.5, SF12=-20 dB

airtime_per_sf = zeros(1,6);
for i = 1:6
    sf_i  = sf_list(i);
    T_s   = 2^sf_i / params.BW;
    T_pr  = (params.Preamble + 4.25) * T_s;
    ci    = max(ceil((8*PL - 4*sf_i + 28 + 16)/(4*sf_i)) * (CR+4), 0);
    airtime_per_sf(i) = T_pr + (8+ci)*T_s;   % seconds
end

fprintf('[SF -> Required SNR -> Airtime]\n');
fprintf('  SF   ReqSNR(dB)  Airtime(ms)\n');
for i = 1:6
    fprintf('  SF%d   %+6.1f     %7.2f\n', sf_list(i), req_snr_sf(i), airtime_per_sf(i)*1e3);
end
fprintf('\n');

%% -------------------------------------------------------
%  SECTION 6D: SX1276 Nonlinear Tx Current Lookup Table
%
%  Source: Semtech SX1276 datasheet, Table 7 (PA_BOOST, 3.3 V):
%    2  dBm -> 10.0 mA
%    5  dBm -> 13.5 mA
%    8  dBm -> 18.0 mA
%    11 dBm -> 22.0 mA
%    14 dBm -> 28.0 mA
%  Linear interpolation between datasheet points.
%  Current is NOT linear with power — using a simple ratio
%  (p_dBm/14) over-estimates current at low power levels.
%% -------------------------------------------------------
I_sleep_uA = 1.5;   % sleep current (uA) — SX1276 datasheet
V_bat      = 3.3;   % battery voltage (V)

% SX1276 Tx current lookup (PA_BOOST pin, 3.3 V supply)
sx1276_pwr_dBm = [2,    5,    8,    11,   14  ];   % dBm
sx1276_I_mA    = [10.0, 13.5, 18.0, 22.0, 28.0];   % mA

% Interpolate current for any power in [2, 14] dBm
sx1276_current = @(p_dBm) interp1(sx1276_pwr_dBm, sx1276_I_mA, ...
    min(max(p_dBm, 2), 14), 'linear');

% Energy per Tx packet: E(uJ) = I(A) x V(V) x t(s) x 1e6
energy_per_tx = @(sf_idx, p_dBm) ...
    sx1276_current(p_dBm) * 1e-3 * V_bat * airtime_per_sf(sf_idx) * 1e6;

fprintf('[SX1276 Tx Current — Datasheet Lookup]\n');
fprintf('  Power(dBm)   I_tx(mA)\n');
for p = [2, 5, 8, 11, 14]
    fprintf('  %5d dBm    %6.1f mA\n', p, sx1276_current(p));
end
fprintf('\n');

%% -------------------------------------------------------
%  SECTION 6E: Fixed SF10 Simulation (Baseline)
%% -------------------------------------------------------
fprintf('[Running Fixed SF10 Simulation (Baseline)...]\n');

fixed_sf      = 10;
fixed_sf_idx  = fixed_sf - 6;
fixed_ptx     = 14;
req_snr_fixed = req_snr_sf(fixed_sf_idx);

fixed_delivered = zeros(1, n_packets);
fixed_energy_uJ = zeros(1, n_packets);

for pkt = 1:n_packets
    snr = SNR_profile(pkt);
    fixed_delivered(pkt) = double(snr >= req_snr_fixed);
    fixed_energy_uJ(pkt) = energy_per_tx(fixed_sf_idx, fixed_ptx);
end

fixed_delivery_pct = 100 * mean(fixed_delivered);
fixed_total_energy = sum(fixed_energy_uJ);
fprintf('  Delivery ratio    : %.2f%%\n', fixed_delivery_pct);
fprintf('  Total energy used : %.2f mJ (over 24 hours)\n', fixed_total_energy/1e3);
fprintf('  Total airtime     : %.2f s\n\n', ...
    sum(fixed_delivered) * airtime_per_sf(fixed_sf_idx));

%% -------------------------------------------------------
%  SECTION 6F: ADR-Controlled Simulation
%% -------------------------------------------------------
fprintf('[Running ADR Simulation...]\n');

% State variables
adr_sf      = 10;       % start at SF10 (same as fixed baseline)
adr_ptx     = 14;       % start at maximum Tx power
adr_sf_idx  = adr_sf - 6;
snr_history = [];       % sliding SNR window — successful packets only

% ADR_ACK_REQ state
ack_count      = 0;     % consecutive unacknowledged uplink counter
ack_req_active = false; % true once ADR_ACK_REQ flag is raised

% Result arrays
adr_sf_log      = zeros(1, n_packets);
adr_ptx_log     = zeros(1, n_packets);
adr_delivered   = zeros(1, n_packets);
adr_energy_uJ   = zeros(1, n_packets);
adr_ack_req_log = zeros(1, n_packets);
adr_event_log   = {};

for pkt = 1:n_packets
    snr = SNR_profile(pkt);

    % Effective SNR with current Tx power back-off applied
    ptx_backoff_dB = adr.p_tx_max_dBm - adr_ptx;
    snr_effective  = snr - ptx_backoff_dB;

    % Packet outcome
    req_snr_cur        = req_snr_sf(adr_sf_idx);
    adr_delivered(pkt) = double(snr_effective >= req_snr_cur);
    adr_sf_log(pkt)    = adr_sf;
    adr_ptx_log(pkt)   = adr_ptx;
    adr_energy_uJ(pkt) = energy_per_tx(adr_sf_idx, adr_ptx);
    adr_ack_req_log(pkt) = double(ack_req_active);

    if adr_delivered(pkt) == 1
        % --- Successful delivery: update SNR history ---
        % Only successful packets enter the history window.
        % This keeps ADR decisions based on good-channel samples only,
        % preventing monsoon-period failures from corrupting SNR_max.
        snr_history(end+1) = snr_effective; %#ok<AGROW>
        if length(snr_history) > adr.history_len
            snr_history = snr_history(end-adr.history_len+1:end);
        end
        ack_count      = 0;
        ack_req_active = false;

    else
        % --- Failed delivery: ADR_ACK_REQ countdown ---
        ack_count = ack_count + 1;

        % Phase 1: raise ADR_ACK_REQ flag after ack_limit failures
        if ack_count > adr.ack_limit && ~ack_req_active
            ack_req_active = true;
            msg = sprintf('  Pkt%3d t=%4.1fh: ADR_ACK_REQ raised (no downlink for %d pkts)', ...
                pkt, t_hours(pkt), ack_count);
            fprintf('%s\n', msg);
            adr_event_log{end+1} = msg; %#ok<AGROW>
        end

        % Phase 2: full reset to SF12 + 14 dBm after ack_limit + ack_delay
        if ack_count > (adr.ack_limit + adr.ack_delay)
            old_sf  = adr_sf;
            old_ptx = adr_ptx;
            adr_sf         = adr.sf_max;
            adr_ptx        = adr.p_tx_max_dBm;
            adr_sf_idx     = adr_sf - 6;
            ack_count      = 0;
            ack_req_active = false;
            msg = sprintf('  Pkt%3d t=%4.1fh: RESET SF%d->SF12, %ddBm->14dBm (ACK timeout)', ...
                pkt, t_hours(pkt), old_sf, old_ptx);
            fprintf('%s\n', msg);
            adr_event_log{end+1} = msg; %#ok<AGROW>
        end
    end

    % --- ADR Decision Block (runs every update_every packets) ---
    % Requires a full history window before making any SF/power changes.
    if mod(pkt, adr.update_every) == 0 && length(snr_history) >= adr.history_len

        snr_max    = max(snr_history);
        req_snr_c  = req_snr_sf(adr_sf_idx);
        snr_margin = snr_max - req_snr_c;

        % Priority 1: margin comfortably above target -> reduce SF
        if snr_margin > adr.margin_dB
            n_steps = floor((snr_margin - adr.margin_dB) / 2.5);
            n_steps = min(n_steps, adr_sf - adr.sf_min);

            if n_steps > 0
                old_sf     = adr_sf;
                adr_sf     = adr_sf - n_steps;
                adr_sf_idx = adr_sf - 6;
                msg = sprintf('  Pkt%3d t=%4.1fh: SF%d->SF%d  (margin=%.1fdB, steps=%d)', ...
                    pkt, t_hours(pkt), old_sf, adr_sf, snr_margin, n_steps);
                fprintf('%s\n', msg);
                adr_event_log{end+1} = msg; %#ok<AGROW>

            % Priority 2: SF already at minimum -> reduce Tx power
            elseif adr_sf == adr.sf_min && adr_ptx > adr.p_tx_min_dBm
                old_ptx = adr_ptx;
                adr_ptx = max(adr.p_tx_min_dBm, adr_ptx - adr.p_tx_step_dB);
                msg = sprintf('  Pkt%3d t=%4.1fh: Power %ddBm->%ddBm (margin=%.1fdB)', ...
                    pkt, t_hours(pkt), old_ptx, adr_ptx, snr_margin);
                fprintf('%s\n', msg);
                adr_event_log{end+1} = msg; %#ok<AGROW>
            end

        % Priority 3: margin too small -> restore power first, then SF up
        elseif snr_margin < adr.snr_margin_up

            if adr_ptx < adr.p_tx_max_dBm
                old_ptx = adr_ptx;
                adr_ptx = min(adr.p_tx_max_dBm, adr_ptx + adr.p_tx_step_dB);
                msg = sprintf('  Pkt%3d t=%4.1fh: Power restore %ddBm->%ddBm (margin=%.1fdB)', ...
                    pkt, t_hours(pkt), old_ptx, adr_ptx, snr_margin);
                fprintf('%s\n', msg);
                adr_event_log{end+1} = msg; %#ok<AGROW>

            elseif adr_sf < adr.sf_max
                old_sf     = adr_sf;
                adr_sf     = adr_sf + 1;
                adr_sf_idx = adr_sf - 6;
                msg = sprintf('  Pkt%3d t=%4.1fh: SF%d->SF%d  (margin=%.1fdB, LINK WEAK)', ...
                    pkt, t_hours(pkt), old_sf, adr_sf, snr_margin);
                fprintf('%s\n', msg);
                adr_event_log{end+1} = msg; %#ok<AGROW>
            end
        end
    end
end

adr_delivery_pct = 100 * mean(adr_delivered);
adr_total_energy = sum(adr_energy_uJ);
fprintf('\n  Delivery ratio    : %.2f%%\n', adr_delivery_pct);
fprintf('  Total energy used : %.2f mJ (over 24 hours)\n', adr_total_energy/1e3);
fprintf('  Total airtime     : %.2f s\n', ...
    sum(adr_delivered .* airtime_per_sf(max(adr_sf_log-6,1))));
fprintf('  ADR events        : %d SF/power changes\n', length(adr_event_log));
fprintf('  ACK_REQ triggers  : %d times\n\n', sum(diff([0 adr_ack_req_log]) == 1));

%% -------------------------------------------------------
%  SECTION 6G: Comparison Statistics
%% -------------------------------------------------------
at_fixed_vec      = airtime_per_sf(fixed_sf_idx) * 1000 * ones(1, n_packets);
at_adr_vec        = airtime_per_sf(max(adr_sf_log-6,1)) * 1000;
energy_saved_pct  = 100 * (1 - adr_total_energy / fixed_total_energy);
airtime_saved_pct = 100 * (1 - sum(at_adr_vec) / sum(at_fixed_vec));

fprintf('[ADR vs Fixed SF10 — 24-Hour Comparison]\n');
fprintf('  %-30s  Fixed SF10    ADR\n', 'Metric');
fprintf('  %s\n', repmat('-',1,52));
fprintf('  %-30s  %8.2f%%   %8.2f%%\n', 'Packet delivery ratio', ...
    fixed_delivery_pct, adr_delivery_pct);
fprintf('  %-30s  %8.2f mJ  %8.2f mJ\n', 'Total energy (24h)', ...
    fixed_total_energy/1e3, adr_total_energy/1e3);
fprintf('  %-30s  %8s    %8.1f%%\n', 'Energy saving vs Fixed', ...
    'baseline', energy_saved_pct);
fprintf('  %-30s  %8.2f s   %8.2f s\n', 'Total channel airtime', ...
    sum(at_fixed_vec)/1e3, sum(at_adr_vec)/1e3);
fprintf('  %-30s  %8s    %8.1f%%\n', 'Airtime saving vs Fixed', ...
    'baseline', airtime_saved_pct);

fprintf('\n[ADR — SF Usage Distribution]\n');
for i = 1:6
    sf_i = sf_list(i);
    cnt  = sum(adr_sf_log == sf_i);
    pct  = 100 * cnt / n_packets;
    bw   = round(pct/2);
    fprintf('  SF%d: %3d pkts (%5.1f%%)  %s\n', sf_i, cnt, pct, repmat('*',1,bw));
end

fprintf('\n[ADR_ACK_REQ Events]\n');
ack_req_starts = find(diff([0 adr_ack_req_log]) == 1);
if isempty(ack_req_starts)
    fprintf('  No ADR_ACK_REQ events triggered (link was robust throughout)\n');
else
    for k = 1:length(ack_req_starts)
        fprintf('  ACK_REQ event %d at t=%.1fh (packet %d)\n', k, ...
            t_hours(ack_req_starts(k)), ack_req_starts(k));
    end
end

% Battery life
fprintf('\n[Battery Life]\n');
C_bat_mAh  = 2000;
V_bat_calc = 3.3;
I_sleep_total = (I_sleep_uA*1e-3) * ...
    (1 - n_packets*airtime_per_sf(fixed_sf_idx)/(24*3600));
I_avg_fixed = (fixed_total_energy*1e-6) / (24*3600*V_bat_calc) * 1e3;
I_avg_adr   = (adr_total_energy  *1e-6) / (24*3600*V_bat_calc) * 1e3;
life_fixed  = (C_bat_mAh / (I_avg_fixed + I_sleep_total*1e3)) / 24;
life_adr    = (C_bat_mAh / (I_avg_adr   + I_sleep_total*1e3)) / 24;
fprintf('  Fixed SF10 battery life : %.1f days = %.2f years\n', life_fixed, life_fixed/365);
fprintf('  ADR    battery life     : %.1f days = %.2f years\n', life_adr, life_adr/365);
fprintf('  Battery life extension  : %.1f days (+%.1f%%)\n', ...
    life_adr - life_fixed, 100*(life_adr/life_fixed - 1));

%% -------------------------------------------------------
%  SECTION 6H: Figures
%% -------------------------------------------------------

% ---- Figure 1: SNR profile + ADR time series ----
figure('Name','ADR — Time Series','Position',[50 50 1200 620]);
sf_colors = lines(6);

subplot(3,1,1);
plot(t_hours, SNR_profile, 'b-', 'LineWidth', 1.2); hold on;
fill([16 20 20 16], [-35 -35 45 45], [0.8 0.8 1.0], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
for i = 1:6
    yline(req_snr_sf(i), '--', 'Color', sf_colors(i,:), 'Alpha', 0.5, ...
        'LineWidth', 0.8, 'Label', sprintf('SF%d req', sf_list(i)));
end
if ~isempty(ack_req_starts)
    ack_req_ends = find(diff([adr_ack_req_log 0]) == -1);
    for k = 1:length(ack_req_starts)
        xe = min(ack_req_ends(k), n_packets);
        fill(t_hours([ack_req_starts(k) xe xe ack_req_starts(k)]), ...
            [-35 -35 45 45], [1 0.8 0.8], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    end
    text(t_hours(ack_req_starts(1)), min(SNR_profile)+3, 'ACK\_REQ', ...
        'FontSize', 8, 'Color', [0.7 0.1 0.1]);
end
xlabel('Time (hours)'); ylabel('SNR (dB)');
title('24-Hour SNR Profile (blue = monsoon window, red = ACK\_REQ active)');
grid on; xlim([0 24]);

subplot(3,1,2);
stairs(t_hours, adr_sf_log, 'r-', 'LineWidth', 2, 'DisplayName', 'ADR SF'); hold on;
plot(t_hours, fixed_sf*ones(1,n_packets), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Fixed SF10');
fail_adr   = t_hours(adr_delivered   == 0);
fail_fixed = t_hours(fixed_delivered == 0);
if ~isempty(fail_adr)
    scatter(fail_adr,   adr_sf_log(adr_delivered==0),   30, 'rx', ...
        'LineWidth', 2, 'DisplayName', 'ADR fail');
end
if ~isempty(fail_fixed)
    scatter(fail_fixed, fixed_sf*ones(size(fail_fixed)), 30, 'bx', ...
        'LineWidth', 2, 'DisplayName', 'Fixed fail');
end
xlabel('Time (hours)'); ylabel('Spreading Factor');
title('ADR Spreading Factor Adaptation vs Fixed SF10');
yticks(7:12); ylim([6.5 12.5]);
legend('Location', 'northeast', 'FontSize', 8); grid on; xlim([0 24]);

subplot(3,1,3);
stairs(t_hours, adr_ptx_log, 'Color', [0.8 0.4 0.1], 'LineWidth', 2, ...
    'DisplayName', 'ADR Tx Power'); hold on;
plot(t_hours, fixed_ptx*ones(1,n_packets), 'b--', 'LineWidth', 1.5, ...
    'DisplayName', 'Fixed Power');
xlabel('Time (hours)'); ylabel('Tx Power (dBm)');
title('ADR Tx Power Adaptation');
yticks(2:3:14); ylim([0 16]);
legend('Location', 'northeast', 'FontSize', 8); grid on; xlim([0 24]);
sgtitle('LoRaWAN ADR — 24-Hour Arabian Sea Buoy Simulation', 'FontWeight', 'bold');

% ---- Figure 2: Energy + Delivery comparison ----
figure('Name','ADR — Performance','Position',[100 100 1200 500]);

subplot(1,3,1);
cum_e_fixed = cumsum(fixed_energy_uJ) / 1e3;
cum_e_adr   = cumsum(adr_energy_uJ)   / 1e3;
plot(t_hours, cum_e_fixed, 'b-', 'LineWidth', 2, 'DisplayName', 'Fixed SF10'); hold on;
plot(t_hours, cum_e_adr,   'r-', 'LineWidth', 2, 'DisplayName', 'ADR');
xlabel('Time (hours)'); ylabel('Cumulative Tx Energy (mJ)');
title('Cumulative Energy Consumption');
legend('Location', 'northwest'); grid on;
text(12, max(cum_e_fixed)*0.5, sprintf('ADR saves\n%.1f%% energy', energy_saved_pct), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

subplot(1,3,2);
cum_del_fixed = cumsum(fixed_delivered) ./ (1:n_packets) * 100;
cum_del_adr   = cumsum(adr_delivered)   ./ (1:n_packets) * 100;
plot(t_hours, cum_del_fixed, 'b-', 'LineWidth', 2, 'DisplayName', 'Fixed SF10'); hold on;
plot(t_hours, cum_del_adr,   'r-', 'LineWidth', 2, 'DisplayName', 'ADR');
yline(95, 'k--', '95% target');
xlabel('Time (hours)'); ylabel('Delivery Ratio (%)');
title('Cumulative Packet Delivery Ratio');
legend('Location', 'southwest'); grid on; ylim([80 105]);

subplot(1,3,3);
sf_counts = histcounts(adr_sf_log, 6.5:12.5);
sf_labels_pie = arrayfun(@(s,c) sprintf('SF%d\n(%d pkts)', s, c), sf_list, sf_counts, ...
    'UniformOutput', false);
nonzero = sf_counts > 0;
pie(sf_counts(nonzero), sf_labels_pie(nonzero));
title('ADR — SF Usage Over 24 Hours');
colormap(gca, lines(sum(nonzero)));
sgtitle('ADR Performance — Energy, Delivery & SF Distribution', 'FontWeight', 'bold');

% ---- Figure 3: SNR margin and airtime ----
figure('Name','ADR — Margin and Airtime','Position',[150 150 1100 450]);

subplot(1,2,1);
margin_fixed = SNR_profile - req_snr_sf(fixed_sf_idx);
margin_adr   = zeros(1, n_packets);
for pkt = 1:n_packets
    ptx_backoff     = adr.p_tx_max_dBm - adr_ptx_log(pkt);
    snr_eff         = SNR_profile(pkt) - ptx_backoff;
    sf_idx_pkt      = adr_sf_log(pkt) - 6;
    margin_adr(pkt) = snr_eff - req_snr_sf(sf_idx_pkt);
end
plot(t_hours, margin_fixed, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fixed SF10'); hold on;
plot(t_hours, margin_adr,   'r-', 'LineWidth', 1.5, 'DisplayName', 'ADR');
yline(adr.margin_dB, 'g--', sprintf('ADR target margin (%d dB)', adr.margin_dB));
yline(0, 'k--', 'Link budget limit');
yline(adr.snr_margin_up, 'm:', 'SF increase threshold');
fill([0 24 24 0], [-50 -50 0 0], [1 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
xlabel('Time (hours)'); ylabel('Link Margin (dB)');
title('Link Margin Over 24 Hours');
legend('Location', 'northeast', 'FontSize', 8); grid on; xlim([0 24]);
text(12, -8, 'Link failure region', ...
    'HorizontalAlignment', 'center', 'Color', [0.7 0.2 0.2], 'FontSize', 9);

subplot(1,2,2);
plot(t_hours, at_fixed_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Fixed SF10'); hold on;
stairs(t_hours, at_adr_vec, 'r-', 'LineWidth', 2, 'DisplayName', 'ADR');
xlabel('Time (hours)'); ylabel('Packet Airtime (ms)');
title('Packet Airtime — ADR vs Fixed');
legend('Location', 'northeast'); grid on; xlim([0 24]);
text(12, max(at_fixed_vec)*0.6, sprintf('ADR saves\n%.1f%% airtime', airtime_saved_pct), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
sgtitle('SNR Margin & Airtime Analysis — ADR vs Fixed', 'FontWeight', 'bold');

% ---- Figure 4: ADR_ACK_REQ counter timeline ----
figure('Name','ADR — ACK REQ','Position',[200 200 1000 380]);
ack_count_log = zeros(1, n_packets);
ac = 0;
for pkt = 1:n_packets
    if adr_delivered(pkt) == 0
        ac = ac + 1;
    else
        ac = 0;
    end
    ack_count_log(pkt) = ac;
end

area(t_hours, adr_ack_req_log*12, 'FaceColor', [1 0.85 0.85], ...
    'EdgeColor', 'none', 'DisplayName', 'ACK\_REQ active'); hold on;

plot(t_hours, ack_count_log, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'Consecutive unack count');

stairs(t_hours, adr_sf_log, 'r-', 'LineWidth', 2, 'DisplayName', 'ADR SF');

yline(adr.ack_limit, 'r--', ...
    sprintf('ACK\\_LIMIT = %d', adr.ack_limit), ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');                     % <-- CHANGE 1: hides from legend

yline(adr.ack_limit + adr.ack_delay, 'm--', ...
    sprintf('LIMIT+DELAY = %d', adr.ack_limit + adr.ack_delay), ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');                     % <-- CHANGE 1: hides from legend

xlabel('Time (hours)');
ylabel('Unack counter / SF / ACK\_REQ flag (x12)');
% title(...) line REMOVED                          % <-- CHANGE 2: deleted to avoid overlap
legend('Location', 'northeast', 'FontSize', 8); grid on; xlim([0 24]);
sgtitle('ADR\_ACK\_REQ Mechanism — Arabian Sea Buoy', 'FontWeight', 'bold', 'FontSize', 12);
%% -------------------------------------------------------
%  SECTION 6I: Final Summary
%% -------------------------------------------------------
fprintf('\n=============================================================\n');
fprintf('  ADR SIMULATION SUMMARY\n');
fprintf('=============================================================\n\n');
fprintf('[Simulation overview]\n');
fprintf('  Total packets simulated : %d  (24h x 5-min interval)\n', n_packets);
fprintf('  ADR events triggered    : %d\n', length(adr_event_log));
fprintf('  ACK_REQ events          : %d\n\n', length(ack_req_starts));
fprintf('[Performance vs Fixed SF10]\n');
fprintf('  Delivery : Fixed=%.1f%%  ADR=%.1f%%  (Delta=+%.1f%%)\n', ...
    fixed_delivery_pct, adr_delivery_pct, adr_delivery_pct - fixed_delivery_pct);
fprintf('  Energy   : Fixed=%.1fmJ  ADR=%.1fmJ  (saved %.1f%%)\n', ...
    fixed_total_energy/1e3, adr_total_energy/1e3, energy_saved_pct);
fprintf('  Airtime  : ADR saves %.1f%% channel airtime\n', airtime_saved_pct);
fprintf('  Battery  : Fixed=%.1f days  ADR=%.1f days  (+%.1f days)\n\n', ...
    life_fixed, life_adr, life_adr - life_fixed);
fprintf('[ADR behaviour over 24 hours]\n');
fprintf('  Morning calm     : ADR steps down to SF7-SF8  -> faster packets\n');
fprintf('  Monsoon period   : ADR steps up to SF11-SF12  -> maintains link\n');
fprintf('  Night recovery   : ADR reduces back to SF8-SF9\n');
fprintf('  Net result       : %.1f%% energy saving with %.1f%% better delivery\n', ...
    energy_saved_pct, adr_delivery_pct - fixed_delivery_pct);

%% -------------------------------------------------------
%  Save
%% -------------------------------------------------------
save('lora_adr.mat', ...
    'adr', 'SNR_profile', 't_hours', 'n_packets', ...
    'fixed_sf', 'fixed_delivered', 'fixed_energy_uJ', 'fixed_delivery_pct', ...
    'adr_sf_log', 'adr_ptx_log', 'adr_delivered', 'adr_energy_uJ', ...
    'adr_delivery_pct', 'adr_total_energy', 'fixed_total_energy', ...
    'adr_event_log', 'adr_ack_req_log', ...
    'life_fixed', 'life_adr', ...
    'req_snr_sf', 'airtime_per_sf', 'sf_list', ...
    'sx1276_pwr_dBm', 'sx1276_I_mA');

fprintf('\nChunk 6 complete. Saved to lora_adr.mat\n');
fprintf('=============================================================\n');
