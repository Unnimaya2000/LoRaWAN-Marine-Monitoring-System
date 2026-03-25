%% =========================================================
%  CHUNK 1 of 5 — Setup + Sensor Data + LoRaWAN Packet
%  CUSAT RUSA: Arabian Sea Buoy — Complete Real-World LoRa Simulation
%  Combines: Code A (LoRaWAN framing, CRC-16) + Code B (sensor encoding)
%  Saves: lora_complete.mat
%% =========================================================
clc; clear; close all;
rng(42);

fprintf('=============================================================\n');
fprintf('  COMPLETE REAL-WORLD LoRa SIMULATION\n');
fprintf('  CUSAT RUSA — Arabian Sea Buoy Monitoring System\n');
fprintf('=============================================================\n\n');
fprintf('--- CHUNK 1: Setup + Sensor + LoRaWAN Packet ---\n\n');

%% -------------------------------------------------------
%  SECTION 1A: LoRa Physical Layer Parameters
%% -------------------------------------------------------
params.SF           = 10;
params.BW           = 125e3;        % 125 kHz
params.CR           = 1;            % coding rate index: 4/(4+CR) = 4/5
params.f_carrier    = 865.1e6;      % India IN865 band (Hz)
params.P_tx_dBm     = 14;           % Tx power (dBm) — SX1276 typical
params.Gtx_dBi      = 2.15;         % Tx dipole antenna gain
params.Grx_dBi      = 2.15;         % Rx dipole antenna gain
params.NF_dB        = 6;            % Receiver noise figure (dB)
params.T_K          = 290;          % Noise temperature (Kelvin)
params.k_B          = 1.38e-23;     % Boltzmann constant
params.Preamble     = 8;            % Preamble symbols
params.Fs           = params.BW * 4; % Oversampling x4

% Derived physical layer values
params.N            = 2^params.SF;              % chips per symbol = 1024
params.T_symbol     = params.N / params.BW;     % symbol duration (s)
params.Rs           = 1 / params.T_symbol;      % symbol rate (Hz)
params.Rb           = params.SF * params.Rs * (4/(4 + params.CR)); % bit rate

fprintf('[LoRa Physical Layer Config]\n');
fprintf('  Spreading Factor : SF%d  (%d chips/symbol)\n', params.SF, params.N);
fprintf('  Bandwidth        : %.0f kHz\n',  params.BW/1e3);
fprintf('  Carrier Freq     : %.1f MHz\n',  params.f_carrier/1e6);
fprintf('  Coding Rate      : 4/%d\n',      4+params.CR);
fprintf('  Symbol Duration  : %.3f ms\n',   params.T_symbol*1e3);
fprintf('  Symbol Rate      : %.2f Hz\n',   params.Rs);
fprintf('  Bit Rate         : %.1f bps\n',  params.Rb);
fprintf('  Antenna Gain Tx  : %.2f dBi\n',  params.Gtx_dBi);
fprintf('  Antenna Gain Rx  : %.2f dBi\n',  params.Grx_dBi);
fprintf('  Noise Figure     : %.1f dB\n',   params.NF_dB);

%% -------------------------------------------------------
%  SECTION 1B: Hardware Impairment Parameters
%  (From Code A — will be applied in Chunk 2)
%% -------------------------------------------------------
hw.ppm_node      = 10;    % Node TCXO accuracy (ppm)
hw.ppm_gateway   = 1;     % Gateway TCXO accuracy (ppm)
hw.CFO_Hz        = params.f_carrier * (hw.ppm_node + hw.ppm_gateway) * 1e-6;
hw.clk_ppm       = 20;    % Sampling clock offset (ppm)
hw.pn_linewidth  = 1e3;   % Oscillator phase noise linewidth (Hz)
hw.rx_filter_loss_dB = 0.5; % SAW filter insertion loss (dB)

fprintf('\n[Hardware Impairment Parameters]\n');
fprintf('  CFO (node+GW)    : ±%.1f Hz  (%.0f+%d ppm @ %.0f MHz)\n', ...
    hw.CFO_Hz, hw.ppm_node, hw.ppm_gateway, params.f_carrier/1e6);
sym_bw = params.BW / params.N;
fprintf('  CFO/Symbol BW    : %.4f  (%.2f%% of bin)\n', ...
    hw.CFO_Hz/sym_bw, 100*hw.CFO_Hz/sym_bw);
fprintf('  Clock offset     : %d ppm\n', hw.clk_ppm);
fprintf('  Phase noise BW   : %.0f Hz\n', hw.pn_linewidth);
pn_penalty = 10*log10(1 + hw.pn_linewidth/sym_bw);
fprintf('  Phase noise SNR penalty (SF%d): %.4f dB\n', params.SF, pn_penalty);
fprintf('  Rx filter loss   : %.1f dB\n', hw.rx_filter_loss_dB);

%% -------------------------------------------------------
%  SECTION 1C: Channel Model Parameters (From Code A)
%% -------------------------------------------------------
channel.K_rice_dB       = 10;   % Rician K-factor — strong sea LOS
channel.K_rice          = 10^(channel.K_rice_dB/10);
channel.sigma_shadow_dB = 6;    % Log-normal shadowing std-dev (dB)
channel.rain_margin_dB  = 3;    % Rain fade margin
channel.sea_state_dB    = 2;    % Sea-state tilt loss
channel.ant_tilt_dB     = 1;    % Buoy antenna tilt
channel.total_margin_dB = channel.rain_margin_dB + ...
                          channel.sea_state_dB   + ...
                          channel.ant_tilt_dB;

% TDL multipath taps (open-sea reflections)
channel.tdl_delays_us  = [0,    3.5,  8.0];   % tap delays (µs)
channel.tdl_gains_dB   = [0,   -12,  -18];    % tap gains (dB)
channel.tdl_gains_lin  = 10.^(channel.tdl_gains_dB/20);

fprintf('\n[Channel Model Parameters — Arabian Sea]\n');
fprintf('  Rician K-factor          : %.1f dB\n', channel.K_rice_dB);
fprintf('  Log-normal shadow σ      : %.1f dB\n', channel.sigma_shadow_dB);
fprintf('  Rain margin              : %.1f dB\n', channel.rain_margin_dB);
fprintf('  Sea-state + tilt margin  : %.1f dB\n', channel.sea_state_dB + channel.ant_tilt_dB);
fprintf('  Total practical margin   : %.1f dB\n', channel.total_margin_dB);
fprintf('  TDL taps: delays [%.1f %.1f %.1f] µs, gains [%d %d %d] dB\n', ...
    channel.tdl_delays_us, channel.tdl_gains_dB);

%% -------------------------------------------------------
%  SECTION 1D: Simulated Arabian Sea Buoy Sensor Readings
%% -------------------------------------------------------
sensor.SST          = 28.4;    % Sea Surface Temperature (°C)
sensor.WaterTemp    = 25.1;    % Water Temperature @ 5m depth (°C)
sensor.Salinity     = 36.2;    % Salinity (PSU)
sensor.Conductivity = 52.8;    % Conductivity (mS/cm)
sensor.Heave        = 1.3;     % Heave (m)
sensor.WaveHeight   = 2.1;     % Significant Wave Height (m)
sensor.Lat          = 10.521;  % Latitude (decimal degrees)
sensor.Lon          = 76.214;  % Longitude (decimal degrees)
sensor.BattVoltage  = 3.82;    % Battery voltage (V)

fprintf('\n[Arabian Sea Buoy Sensor Readings]\n');
fprintf('  Sea Surface Temp : %.1f °C\n',   sensor.SST);
fprintf('  Water Temp @5m   : %.1f °C\n',   sensor.WaterTemp);
fprintf('  Salinity         : %.1f PSU\n',  sensor.Salinity);
fprintf('  Conductivity     : %.1f mS/cm\n',sensor.Conductivity);
fprintf('  Heave            : %.1f m\n',    sensor.Heave);
fprintf('  Wave Height      : %.1f m\n',    sensor.WaveHeight);
fprintf('  Position         : %.3f°N, %.3f°E\n', sensor.Lat, sensor.Lon);
fprintf('  Battery          : %.2f V\n',    sensor.BattVoltage);

%% -------------------------------------------------------
%  SECTION 1E: Encode Sensor Values into Bytes
%% -------------------------------------------------------
% Each field encoded to fit in uint8 (0-255) with defined resolution
SST_enc   = uint8(round(sensor.SST   * 2));    % 0.5°C resolution
WTemp_enc = uint8(round(sensor.WaterTemp * 2)); % 0.5°C resolution
Sal_enc   = uint8(round(sensor.Salinity));      % 1 PSU resolution
Cond_enc  = uint8(round(sensor.Conductivity));  % 1 mS/cm resolution
Heave_enc = uint8(round(sensor.Heave * 10));    % 0.1m resolution
Wave_enc  = uint8(round(sensor.WaveHeight * 10));% 0.1m resolution
Batt_enc  = uint8(round(sensor.BattVoltage * 20)); % 0.05V resolution

payload_bytes = [SST_enc, WTemp_enc, Sal_enc, Cond_enc, ...
                 Heave_enc, Wave_enc, Batt_enc];

fprintf('\n[Sensor Payload Encoding  (7 bytes)]\n');
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 0.5°C)\n',  'SST:',   sensor.SST,          SST_enc);
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 0.5°C)\n',  'WTemp:', sensor.WaterTemp,     WTemp_enc);
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 1 PSU)\n',  'Salin:',  sensor.Salinity,    Sal_enc);
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 1 mS/cm)\n','Cond:',   sensor.Conductivity, Cond_enc);
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 0.1m)\n',   'Heave:',  sensor.Heave,        Heave_enc);
fprintf('  %-16s raw=%.1f  → enc=%d  (res: 0.1m)\n',   'Wave:',   sensor.WaveHeight,   Wave_enc);
fprintf('  %-16s raw=%.2f → enc=%d  (res: 0.05V)\n',  'Batt:',   sensor.BattVoltage,  Batt_enc);

%% -------------------------------------------------------
%  SECTION 1F: CRC-16/IBM Checksum (From Code A)
%% -------------------------------------------------------
crc_val = crc16_calc(payload_bytes);
fprintf('\n[CRC-16/IBM on payload]\n');
fprintf('  Payload hex : '); fprintf('%02X ', payload_bytes); fprintf('\n');
fprintf('  CRC-16      : 0x%04X\n', crc_val);

%% -------------------------------------------------------
%  SECTION 1G: Build Full LoRaWAN-like Frame (From Code A)
%  MHDR(1) + DevAddr(4) + FCtrl(1) + FCnt(2) + FPort(1) + Payload(N) + MIC(4)
%% -------------------------------------------------------
MHDR     = uint8(0x40);                              % Unconfirmed Data Up
DevAddr  = [uint8(0xAB), uint8(0xCD), uint8(0x01), uint8(0x23)]; % 4-byte address
FCtrl    = uint8(0x00);                              % ADR=0, FOptsLen=0
FCnt     = [uint8(0x00), uint8(0x01)];               % Frame counter = 1
FPort    = uint8(0x02);                              % App port 2 (sensor data)
% MIC placeholder (4 bytes) — real system uses AES-CMAC with NwkSKey
crc_dbl  = double(crc_val);
MIC      = [uint8(bitand(crc_dbl, 255)), uint8(bitshift(crc_dbl, -8)), ...
            uint8(0xBE), uint8(0xEF)];

% Assemble full frame
packet   = [MHDR, DevAddr, FCtrl, FCnt, FPort, payload_bytes, MIC];

fprintf('\n[LoRaWAN Uplink Frame]\n');
fprintf('  MHDR     (1B): 0x%02X  — Unconfirmed Data Up\n', MHDR);
fprintf('  DevAddr  (4B): '); fprintf('%02X ', DevAddr); fprintf('\n');
fprintf('  FCtrl    (1B): 0x%02X\n', FCtrl);
fprintf('  FCnt     (2B): 0x%02X%02X  — Frame #1\n', FCnt(1), FCnt(2));
fprintf('  FPort    (1B): 0x%02X  — App port 2\n', FPort);
fprintf('  Payload  (%dB): ', length(payload_bytes)); fprintf('%02X ', payload_bytes); fprintf('\n');
fprintf('  MIC      (4B): '); fprintf('%02X ', MIC); fprintf('\n');
fprintf('  ─────────────────────────────────────────\n');
fprintf('  Total frame  : %d bytes  (%d B overhead + %d B payload)\n', ...
    length(packet), length(packet)-length(payload_bytes), length(payload_bytes));

%% -------------------------------------------------------
%  SECTION 1H: Convert Packet to Bit Stream
%% -------------------------------------------------------
bit_matrix = de2bi(packet, 8, 'left-msb');
bit_stream = reshape(bit_matrix', 1, []);

fprintf('\n[Bit Stream]\n');
fprintf('  Packet bytes : %d\n',   length(packet));
fprintf('  Total bits   : %d\n',   length(bit_stream));
fprintf('  First 32 bits: ');
fprintf('%d', bit_stream(1:32)); fprintf('\n');

%% -------------------------------------------------------
%  SECTION 1I: Compute LoRa Airtime for this packet
%% -------------------------------------------------------
PL       = length(payload_bytes);
SF       = params.SF;
BW       = params.BW;
CR       = params.CR;
T_sym    = params.T_symbol;

T_preamble = (params.Preamble + 4.25) * T_sym;
n_payload  = 8 + max(ceil((8*PL - 4*SF + 28 + 16) / (4*SF)) * (CR+4), 0);
T_payload  = n_payload * T_sym;
airtime_ms = (T_preamble + T_payload) * 1e3;

fprintf('\n[Packet Airtime]\n');
fprintf('  Payload      : %d bytes\n',   PL);
fprintf('  Preamble sym : %.2f\n',       params.Preamble + 4.25);
fprintf('  Payload sym  : %d\n',         n_payload);
fprintf('  Symbol dur   : %.3f ms\n',    T_sym*1e3);
fprintf('  Airtime      : %.2f ms\n',    airtime_ms);
fprintf('  Duty cycle   : %.4f%%  (5-min interval)\n', airtime_ms/(300e3)*100);

%% -------------------------------------------------------
%  SECTION 1J: Summary Plot — Packet Structure
%% -------------------------------------------------------
figure('Name','Chunk 1 — Packet Structure','Position',[50 50 900 350]);

% Frame field sizes
fields     = {'MHDR','DevAddr','FCtrl','FCnt','FPort','Payload','MIC'};
field_sizes= [1, 4, 1, 2, 1, length(payload_bytes), 4];
field_colors = [0.2 0.5 0.9; 0.9 0.4 0.2; 0.3 0.7 0.3; ...
                0.8 0.2 0.8; 0.9 0.7 0.1; 0.2 0.7 0.7; 0.7 0.3 0.3];

% Draw frame as stacked bars
subplot(1,2,1);
cum = 0;
for i = 1:length(fields)
    h = barh(1, field_sizes(i), 'stacked', 'FaceColor', field_colors(i,:), 'EdgeColor','w', 'LineWidth',1.5);
    hold on;
    text(cum + field_sizes(i)/2, 1, sprintf('%s\n%dB', fields{i}, field_sizes(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'FontWeight','bold');
    cum = cum + field_sizes(i);
end
set(gca,'YTick',[],'XLabel',[]); xlim([0 cum+0.5]);
xlabel('Bytes'); title('LoRaWAN Frame Structure');
grid on;

% Encoding resolution bar chart
subplot(1,2,2);
enc_vals   = double(payload_bytes);
bar_cols   = field_colors(1:7,:);
b = bar(1:7, enc_vals, 0.6, 'FaceColor','flat');
for i=1:7; b.CData(i,:) = bar_cols(i,:); end
names = {'SST','WTemp','Salin.','Cond.','Heave','Wave','Batt'};
set(gca,'XTick',1:7,'XTickLabel',names,'FontSize',9); xtickangle(30);
ylabel('Encoded value (uint8)');
title('Payload: Encoded Sensor Values');
for i=1:7
    text(i, enc_vals(i)+2, sprintf('%d',enc_vals(i)), ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
grid on;

sgtitle('LoRaWAN Frame Build — CUSAT RUSA Buoy','FontWeight','bold');

%% -------------------------------------------------------
%  Save workspace
%% -------------------------------------------------------
% Assign standalone convenience variables (also stored inside params struct)
N  = params.N;
SF = params.SF;
BW = params.BW;
CR = params.CR;

save('lora_complete.mat', 'params', 'hw', 'channel', 'sensor', ...
     'packet', 'payload_bytes', 'bit_stream', 'crc_val', ...
     'airtime_ms', 'n_payload', 'PL', 'N', 'SF', 'BW', 'CR');

fprintf('\n✓ Chunk 1 complete. Saved to lora_complete.mat\n');
fprintf('  → Run Chunk2_Modulation.m next\n\n');

%% -------------------------------------------------------
%  LOCAL HELPER: CRC-16/IBM
%% -------------------------------------------------------
function crc = crc16_calc(data)
    crc  = uint16(0);
    poly = uint16(0x8005);
    for i = 1:length(data)
        crc = bitxor(crc, bitshift(uint16(data(i)), 8));
        for bit = 1:8
            if bitand(crc, uint16(0x8000))
                crc = bitxor(bitshift(crc, 1), poly);
            else
                crc = bitshift(crc, 1);
            end
        end
    end
end
