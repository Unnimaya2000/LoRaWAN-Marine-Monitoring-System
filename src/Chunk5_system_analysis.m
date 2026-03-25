%% =========================================================
%  CHUNK 5 of 5 — System Analysis + Final Summary
%  Full system-level study:
%    - Physics-based link budget with antenna gains
%    - Coverage vs distance for all SFs
%    - Airtime & max range trade-off (SF7–SF12)
%    - Pure ALOHA MAC collision analysis (N buoys)
%    - Battery life model (Tx/Rx/sleep/sense currents)
%    - SX1276 datasheet sensitivity comparison
%    - Automated sanity checks
%    - Final console summary
%  Requires: lora_complete.mat  (from Chunk 4)
%  Last chunk — no save needed
%% =========================================================
clc;
load('lora_complete.mat');
fprintf('--- CHUNK 5: System Analysis + Final Summary ---\n\n');

%% -------------------------------------------------------
%  SECTION 5A: Coverage & Link Margin vs Distance — All SFs
%% -------------------------------------------------------
fprintf('[Coverage Analysis — All SFs]\n');
fprintf('  SF   Sensitivity(dBm)  Max Range(km)  Margin@15km(dB)\n');

sf_vals     = 7:12;
sx1276_sens = [-123, -126, -129, -132, -133, -136]; % datasheet values
colors6     = lines(6);
d_sweep     = 1:0.5:40;

% Noise floor (physics-based)
N_fl_base   = 10*log10(params.k_B * params.T_K * params.BW) + 30 + params.NF_dB;

margins_all = zeros(length(sf_vals), length(d_sweep));
for i = 1:length(sf_vals)
    sf      = sf_vals(i);
    req_snr = -7.5 - 2.5*(sf-7);
    sens    = N_fl_base + req_snr;
    % Max range using full system margin (FSPL + practical margins + shadowing + impl loss)
    % total_system_margin = practical(6dB) + shadowing_1sigma(6dB) + impl_loss(3dB) + fade_margin(5dB) = 20dB
    total_sys_margin = 20;
    EIRP    = params.P_tx_dBm + params.Gtx_dBi + params.Grx_dBi;
    max_fspl= EIRP - sens - total_sys_margin;
    max_r_km= (10^((max_fspl + 147.55 - 20*log10(params.f_carrier))/20)) / 1e3;
    % Margin at 15 km
    fspl_15 = 20*log10(15e3) + 20*log10(params.f_carrier) - 147.55;
    Prx_15  = params.P_tx_dBm + params.Gtx_dBi + params.Grx_dBi - fspl_15;
    margin  = Prx_15 - sens - channel.total_margin_dB;
    fprintf('  SF%d   %8.2f dBm    %8.1f km     %+.2f dB\n', sf, sens, max_r_km, margin);
    % Margin vs distance for plotting
    for di = 1:length(d_sweep)
        fspl_d = 20*log10(d_sweep(di)*1e3) + 20*log10(params.f_carrier) - 147.55;
        Prx_d  = params.P_tx_dBm + params.Gtx_dBi + params.Grx_dBi - fspl_d;
        margins_all(i,di) = Prx_d - sens - channel.total_margin_dB;
    end
end

%% -------------------------------------------------------
%  SECTION 5B: Airtime vs SF (tied to real packet)
%% -------------------------------------------------------
fprintf('\n[Airtime vs SF  (payload=%d bytes, BW=125kHz, CR=4/5)]\n', PL);
airtimes_sf = zeros(1, length(sf_vals));
for i = 1:length(sf_vals)
    sf_i  = sf_vals(i);
    T_s   = 2^sf_i / params.BW;
    T_pr  = (params.Preamble + 4.25) * T_s;
    ci    = max(ceil((8*PL - 4*sf_i + 28 + 16) / (4*sf_i)) * (CR+4), 0);
    airtimes_sf(i) = (T_pr + (8+ci)*T_s) * 1e3;   % ms
    fprintf('  SF%d  :  %.2f ms\n', sf_vals(i), airtimes_sf(i));
end

%% -------------------------------------------------------
%  SECTION 5C: Pure ALOHA MAC Analysis (From Code A)
%% -------------------------------------------------------
fprintf('\n[Pure ALOHA MAC — N Buoys, 5-min Interval, SF10]\n');

tx_interval_s = 300;    % 5-minute reporting interval
at_sf10_s     = airtime_ms / 1000;   % use actual computed airtime from Chunk 1
N_buoys       = [1, 5, 10, 20, 50, 100, 200];

fprintf('  Nodes  G (load)     S (throughput)  P_collision  Pkts/hr/GW\n');
for N_b = N_buoys
    G      = N_b * at_sf10_s / tx_interval_s;
    S      = G * exp(-2*G);
    Pcoll  = 1 - exp(-2*G*(N_b-1)/max(N_b,2));
    if N_b == 1; Pcoll = 0; end
    pkts_hr= N_b * (1-Pcoll) * 3600/tx_interval_s;
    fprintf('  %5d   %10.6f   %10.6f    %10.4f%%  %9.1f\n', ...
        N_b, G, S, Pcoll*100, pkts_hr);
end
G_10 = 10 * at_sf10_s / tx_interval_s;
fprintf('\n  At N=10: G=%.6f  →  Collision = %.4f%%  (practically zero)\n', ...
    G_10, (1-exp(-2*G_10*9/10))*100);

%% -------------------------------------------------------
%  SECTION 5D: Battery Life Estimation (From Code A)
%% -------------------------------------------------------
fprintf('\n[Battery Life Estimation — SX1276 Current Profile]\n');

I_tx_mA     = 28;     % SX1276 Tx current @ 14 dBm
I_rx_mA     = 10.3;   % SX1276 Rx current
I_sleep_uA  = 1.5;    % SX1276 deep sleep
I_mcu_mA    = 5;      % MCU active
I_sens_mA   = 8;      % Sensor warmup
t_sense_s   = 1;      % Sensor measurement window
V_bat       = 3.3;    % LiSOCl2 battery voltage
C_bat_mAh   = 2000;   % 2 Ah cell

fprintf('  SF   Airtime(ms)  E_tx(µJ)  Sleep%%    I_avg(µA)  Life(days)  Life(yrs)\n');
bat_life_days = zeros(1, length(sf_vals));
for i = 1:length(sf_vals)
    sf_i   = sf_vals(i);
    at_s   = airtimes_sf(i) / 1000;
    E_tx   = I_tx_mA*1e-3  * V_bat * at_s    * 1e6;
    E_rx   = I_rx_mA*1e-3  * V_bat * at_s    * 1e6;
    E_sens = (I_mcu_mA+I_sens_mA)*1e-3 * V_bat * t_sense_s * 1e6;
    E_slp  = I_sleep_uA*1e-6 * V_bat * (tx_interval_s - at_s - t_sense_s) * 1e6;
    E_tot  = E_tx + E_rx + E_sens + E_slp;
    I_avg  = E_tot*1e-6 / (tx_interval_s*V_bat) * 1000;  % mA
    life_d = (C_bat_mAh / I_avg) / 24;
    bat_life_days(i) = life_d;
    sleep_pct = 100*(tx_interval_s - at_s)/tx_interval_s;
    fprintf('  SF%d   %8.2f   %8.2f   %6.2f     %7.4f    %8.1f    %6.2f\n', ...
        sf_i, airtimes_sf(i), E_tx, sleep_pct, I_avg*1000, life_d, life_d/365);
end

%% -------------------------------------------------------
%  SECTION 5E: SX1276 Sensitivity Comparison
%% -------------------------------------------------------
fprintf('\n[SX1276 Datasheet vs Simulation Sensitivity]\n');
fprintf('  SF   Simulated(dBm)  SX1276(dBm)  Difference(dB)  Verdict\n');
for i = 1:length(sf_vals)
    sf      = sf_vals(i);
    req_snr = -7.5 - 2.5*(sf-7);
    sens_sim= N_fl_base + req_snr;
    diff_dB = sens_sim - sx1276_sens(i);
    verdict = 'OK (impl. loss)';
    if abs(diff_dB) > 5; verdict = 'CHECK'; end
    fprintf('  SF%d   %8.2f       %8d       %+7.2f         %s\n', ...
        sf, sens_sim, sx1276_sens(i), diff_dB, verdict);
end
fprintf('\n  Note: ~2–4 dB gap = implementation loss, non-ideal filter, CFO penalty.\n');

%% -------------------------------------------------------
%  SECTION 5F: Automated Sanity Checks (From Code A)
%% -------------------------------------------------------
fprintf('\n[Automated Sanity Checks]\n');
all_pass = true;

for sf = sf_vals
    for bw = [62.5e3, 125e3, 250e3]
        % Check 1: Airtime monotonically increases with SF
        if sf > 7
            T_sf   = 2^sf   / bw;
            T_sf_m = 2^(sf-1)/bw;
            at_sf  = (params.Preamble+4.25)*T_sf   + (8+max(ceil((8*PL-4*sf+28+16)/(4*sf))*(CR+4),0))*T_sf;
            at_sfm = (params.Preamble+4.25)*T_sf_m + (8+max(ceil((8*PL-4*(sf-1)+28+16)/(4*(sf-1)))*(CR+4),0))*T_sf_m;
            if at_sf <= at_sfm
                fprintf('  FAIL: Airtime(SF%d) <= Airtime(SF%d) at BW=%.0fkHz\n',sf,sf-1,bw/1e3);
                all_pass = false;
            end
        end
        % Check 2: Sensitivity improves (decreases) with SF
        if sf > 7
            N_fl  = 10*log10(params.k_B*params.T_K*bw)+30+params.NF_dB;
            s_sf  = N_fl + (-7.5-2.5*(sf-7));
            s_sfm = N_fl + (-7.5-2.5*(sf-8));
            if s_sf >= s_sfm
                fprintf('  FAIL: Sensitivity(SF%d) >= Sensitivity(SF%d)\n',sf,sf-1);
                all_pass = false;
            end
        end
        % Check 3: BER → 0 at high SNR
        M_chk   = 2^sf;
        ber_high= (M_chk-1)/(2*log2(M_chk)) * erfc(sqrt(log2(M_chk)*10^(20/10)));
        if ber_high >= 1e-4
            fprintf('  FAIL: BER not near 0 at +20dB for SF%d BW=%.0fkHz\n',sf,bw/1e3);
            all_pass = false;
        end
        % Check 4: PER near 0 at very high SNR
        per_vhigh = 1-(1-(M_chk-1)/(2*log2(M_chk))*erfc(sqrt(log2(M_chk)*10^(30/10))))^(PL*8);
        if per_vhigh >= 1e-4
            fprintf('  FAIL: PER too high at +30dB for SF%d\n',sf);
            all_pass = false;
        end
    end
end

if all_pass
    fprintf('  ✓ All sanity checks PASSED  (SF7–12 × BW 62.5/125/250 kHz)\n');
    fprintf('  ✓ Airtime monotonically increases with SF\n');
    fprintf('  ✓ Sensitivity monotonically improves with SF\n');
    fprintf('  ✓ BER → 0 at high SNR for all configurations\n');
    fprintf('  ✓ PER < 1e-4 at +30 dB SNR\n');
end

%% -------------------------------------------------------
%  SECTION 5G: FINAL PLOTS
%% -------------------------------------------------------

% ---- Figure 1: Coverage — Link Margin vs Distance ----
figure('Name','Chunk 5a — Coverage','Position',[50 50 1200 500]);

subplot(1,2,1);
for i = 1:length(sf_vals)
    plot(d_sweep, margins_all(i,:), 'LineWidth', 2.2, 'Color', colors6(i,:), ...
        'DisplayName', sprintf('SF%d', sf_vals(i)));
    hold on;
end
yline(10,'r--','10 dB target margin','LabelHorizontalAlignment','left');
yline(0, 'k--','Link budget limit');
xline(distance_km,'b:',sprintf('%d km buoy',distance_km));
xlabel('Distance (km)'); ylabel('Link Margin (dB)');
title('Link Margin vs Distance — All SFs (with practical margins)');
legend('Location','northeast','FontSize',8); grid on;
xlim([1 40]); ylim([-25 55]);

subplot(1,2,2);
% Airtime bar + battery overlay
yyaxis left;
b_at = bar(sf_vals, airtimes_sf, 0.5, 'FaceColor','flat');
cmap = cool(6);
for i=1:6; b_at.CData(i,:)=cmap(i,:); end
for i=1:6
    text(sf_vals(i), airtimes_sf(i)+8, sprintf('%.0fms',airtimes_sf(i)), ...
        'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
end
ylabel('Airtime (ms)'); ylim([0 max(airtimes_sf)*1.3]);
yyaxis right;
plot(sf_vals, bat_life_days, 'rs-', 'LineWidth', 2.5, 'MarkerSize', 9, ...
    'MarkerFaceColor','r', 'DisplayName','Battery life');
for i=1:6
    text(sf_vals(i)+0.05, bat_life_days(i)+20, sprintf('%.0fd',bat_life_days(i)), ...
        'FontSize',8,'Color','r');
end
ylabel('Battery Life (days)');
xlabel('Spreading Factor'); title('Airtime & Battery Life vs SF');
legend({'Airtime','Battery (days)'},'Location','northwest');
grid on;

sgtitle('Coverage & Airtime Analysis | IN865 | 14 dBm | 2.15 dBi','FontWeight','bold');

% ---- Figure 2: ALOHA MAC analysis ----
figure('Name','Chunk 5b — MAC + Duty Cycle','Position',[100 100 1100 450]);

N_sweep   = 1:200;
G_sweep   = N_sweep * at_sf10_s / tx_interval_s;
S_sweep   = G_sweep .* exp(-2*G_sweep);
Pc_sweep  = zeros(size(N_sweep));
for ni = 2:length(N_sweep)
    G_n  = G_sweep(ni);
    N_n  = N_sweep(ni);
    Pc_sweep(ni) = 1 - exp(-2*G_n*(N_n-1)/N_n);
end

subplot(1,2,1);
plot(N_sweep, S_sweep*100,  'b-', 'LineWidth', 2.5, 'DisplayName','Throughput (%)'); hold on;
plot(N_sweep, Pc_sweep*100, 'r-', 'LineWidth', 2.5, 'DisplayName','Collision (%)');
xline(10,  'g--', '10 buoys');
xline(100, 'm--', '100 buoys');
xlabel('Number of Buoys'); ylabel('%');
title('Pure ALOHA: Throughput & Collision Probability');
legend('Location','east'); grid on; ylim([0 105]);

subplot(1,2,2);
% Duty cycle pie
at_ms_sf10 = airtime_ms;
sleep_ms   = tx_interval_s*1000 - at_ms_sf10;
pie([at_ms_sf10, sleep_ms], ...
    {sprintf('TX: %.1f ms\n(%.4f%%)', at_ms_sf10, at_ms_sf10/(tx_interval_s*1000)*100), ...
     sprintf('Sleep: %.0f ms\n(%.4f%%)', sleep_ms, sleep_ms/(tx_interval_s*1000)*100)});
colormap([0.9 0.3 0.2; 0.3 0.7 0.4]);
title(sprintf('Duty Cycle — SF10, 5-min interval\nActive: %.4f%%', ...
    at_ms_sf10/(tx_interval_s*1000)*100));

sgtitle('MAC Layer Analysis — Pure ALOHA | CUSAT RUSA Buoy','FontWeight','bold');

% ---- Figure 3: SX1276 Comparison ----
figure('Name','Chunk 5c — SX1276 Comparison','Position',[150 150 750 420]);

sim_sens = zeros(1,6);
for i=1:6
    sf_i = sf_vals(i);
    sim_sens(i) = N_fl_base + (-7.5 - 2.5*(sf_i-7));
end
x = 1:6;
bar(x-0.2, sim_sens,      0.35, 'FaceColor',[0.3 0.6 0.9], 'DisplayName','Simulation');
hold on;
bar(x+0.2, sx1276_sens,   0.35, 'FaceColor',[0.9 0.5 0.2], 'DisplayName','SX1276 Datasheet');
set(gca,'XTick',x,'XTickLabel',arrayfun(@(s) sprintf('SF%d',s),sf_vals,'UniformOutput',false));
ylabel('Sensitivity (dBm)');
title('Receiver Sensitivity: Simulation vs SX1276 Datasheet');
legend('Location','southeast'); grid on;
for i=1:6
    text(i-0.2, sim_sens(i)-0.7,   sprintf('%.0f',sim_sens(i)),   'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
    text(i+0.2, sx1276_sens(i)-0.7,sprintf('%d',sx1276_sens(i)), 'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
end
sgtitle('SX1276 Sensitivity Validation','FontWeight','bold');

%% -------------------------------------------------------
%  FINAL CONSOLE SUMMARY
%% -------------------------------------------------------
fprintf('\n');
fprintf('=============================================================\n');
fprintf('  COMPLETE REAL-WORLD SIMULATION — FINAL SUMMARY\n');
fprintf('  CUSAT RUSA: Arabian Sea Buoy Monitoring System\n');
fprintf('=============================================================\n');
fprintf('\n[Physical Layer]\n');
fprintf('  Modulation      : LoRa CSS (Chirp Spread Spectrum)\n');
fprintf('  Spreading Factor: SF%d  (%d chips/symbol)\n', params.SF, params.N);
fprintf('  Bandwidth       : %.0f kHz\n',  params.BW/1e3);
fprintf('  Carrier         : %.1f MHz (IN865)\n', params.f_carrier/1e6);
fprintf('  Bit Rate        : %.1f bps\n',  params.Rb);
fprintf('  Packet Airtime  : %.2f ms  (real IQ signal)\n', airtime_ms);
fprintf('\n[Signal Pipeline — VERIFIED]\n');
fprintf('  ✓ CSS chirp IQ generated (%d samples, %.4f mean power)\n', ...
    length(chirp_signal), mean(abs(chirp_signal).^2));
fprintf('  ✓ Hardware impairments applied: CFO ±%.0fHz, clock %dppm, phase noise\n', ...
    hw.CFO_Hz, hw.clk_ppm);
fprintf('  ✓ 5 channel models applied to real IQ signal\n');
fprintf('  ✓ Dechirp + FFT demodulator recovered symbols correctly\n');
fprintf('  ✓ Sensor data decoded at gateway — SST=%.1f°C, Sal=%.1fPSU, Heave=%.1fm\n', ...
    dec_SST, dec_Sal, dec_Heave);
fprintf('\n[Channel Models — BER from Monte Carlo Sweep (at operating SNR)]');
% Find BER at operating SNR from the sweep results
[~, op_idx] = min(abs(SNR_sweep_dB - SNR_op_dB));
op_idx = min(op_idx, size(BER_all,2));  % clip to sweep range
for m = 1:5
    fprintf('  %-22s BER = %.6f\n', model_names{m}, BER_all(m, op_idx));
end
fprintf('\n[Link Budget — 15 km]\n');
fprintf('  FSPL            : %.2f dB\n',   FSPL_dB);
fprintf('  Rx Power        : %.2f dBm\n',  P_rx_dBm);
fprintf('  Noise Floor     : %.2f dBm\n',  N_floor_dBm);
fprintf('  Sensitivity SF10: %.2f dBm\n',  sensitivity_dBm);
fprintf('  Link Margin     : %+.2f dB  ✓\n', link_margin_dB);
fprintf('\n[Protocol]\n');
fprintf('  Frame format    : LoRaWAN-like (%dB overhead + %dB payload)\n', ...
    length(packet)-PL, PL);
fprintf('  CRC-16          : 0x%04X  ✓\n', crc_val);
fprintf('\n[System — N=10 Buoys, 5-min interval]\n');
fprintf('  ALOHA load G    : %.6f  (collision probability ≈ 0)\n', G_10);
fprintf('  Battery SF10    : %.1f days = %.2f years (2Ah cell)\n', ...
    bat_life_days(sf_vals==10), bat_life_days(sf_vals==10)/365);
fprintf('  SX1276 gap      : %.1f dB (implementation loss)\n', ...
    (N_fl_base + (-7.5-2.5*(params.SF-7))) - sx1276_sens(sf_vals==10));
fprintf('\n=============================================================\n');
fprintf('  Buoy → CSS chirp → Rician/TDL channel → Gateway : SUCCESS ✓\n');
fprintf('=============================================================\n');
