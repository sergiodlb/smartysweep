%% test HEMT by sweeping gate and measuring Vout with fixed Vdd
transistor_fields.columns = {'V_g (V)', 'V_{amp} (V)'}; % cannot read Vdd since SeekatUNO lacks read operation
transistor_fields.channels = {'Vg', 'K2000.V'};
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(1, 'transistor_test_sweep_testconfigA_Vdd0.2', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;

%% transistor characteristic looked strange; try floating both gates
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(2, 'transistor_test_sweep_testconfigA_Vdd0.2_float_gates', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;
% this works; Vg* = -0.177

%% float only GBG
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(3, 'transistor_test_sweep_testconfigA_Vdd0.2_float_GBG', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;

%% float only BG
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(4, 'transistor_test_sweep_testconfigA_Vdd0.2_float_BG', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;
% this works; Vg* = -0.177

%% float only TG
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(5, 'transistor_test_sweep_testconfigA_Vdd0.2_float_TG', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;

%% BG connected to K2400 through 10 MOhm sourcing zero volts
smset('Vdd', 0.2, 1); % Vdd
smset('Vg', 0, 1);    % Vg
pause(3);
Vg_col = 1;
DC_voltage_sweep(6, 'transistor_test_sweep_testconfigA_Vdd0.2_Vbg0', 0, -0.3, 0.0, 301, Vg_col, transistor_fields, {2});
transistorOFF;
% unless I can "lift" the ground on the BG voltage (to prevent leakage), it
% seems that I will have to float the BG in order to get meaningful
% measurements with the device connected this way

%% measure the standard capacitor with BG grounded
configSTD.Vex_amplitude_channel = 'ZI.OUT1AMP1';
configSTD.Vex_range_channel     = 'ZI.OUT1RANGE';
configSTD.Vstd_amplitude_channel= 'ZI.OUT2AMP2';
configSTD.Vstd_range_channel    = 'ZI.OUT2RANGE';
configSTD.Vstd_phase_channel    = 'ZI.2phase';
configSTD.time_constant_channel = 'ZI.1timeConstant';
configSTD.frequency_channel     = 'ZI.1freq';
configSTD.Rchip                 = 100e6;
configSTD.channels = {'ZI.1X', 'ZI.1Y'};
configSTD.Xcol = 1;
configSTD.Ycol = 2;
transistorON;

[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 150, 21, configSTD, 'Vex', 1e-3, 'Vstd_range', 1e-3);
save('freq_sweep_after_006_BG_ground.mat', 'freq', 'Vc0Vex', 'Vr0Vex');

%% repeat measurement with BG floating
[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 150, 21, configSTD, 'Vex', 1e-3, 'Vstd_range', 1e-3);
save('freq_sweep_after_006_BG_float.mat', 'freq', 'Vc0Vex', 'Vr0Vex');

%% repeat measurement with BG floating and larger Vex
[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 150, 21, configSTD, 'Vex', 9e-3, 'Vstd_range', 10e-3);
save('freq_sweep_after_006_BG_float_Vex9mV.mat', 'freq', 'Vc0Vex', 'Vr0Vex');

%% repeat measurement with BG floating and even larger Vex
[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 150, 21, configSTD, 'Vex', 20e-3, 'Vstd_range', 100e-3);
save('freq_sweep_after_006_BG_float_Vex20mV.mat', 'freq', 'Vc0Vex', 'Vr0Vex');
% why are these measurements scaling with Vex??

%% repeat measurement with BG floating and even larger Vex
[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(40, 200, 21, configSTD, 'Vex', 10e-3, 'Vstd_range', 10e-3);
save('freq_sweep_after_006_BG_float_Vex10mV.mat', 'freq', 'Vc0Vex', 'Vr0Vex');
% *actually they change mostly at low excitation voltage

%% plot together
files = {'freq_sweep_after_006_BG_float.mat', 'freq_sweep_after_006_BG_float_Vex9mV.mat', 'freq_sweep_after_006_BG_float_Vex20mV.mat'};
colors = 'kbr';
f_cutoff = 40;
figure();
for n = 1:3
    load(files{n});
    plot(freq, Vr0Vex, [colors(n), 'o-']);
    hold on;
    plot(freq, Vc0Vex, [colors(n), '--']);
    Cstd = -mean(diff(Vr0Vex(freq>f_cutoff))/diff(freq(freq>f_cutoff)))/(2*pi*100e6)*1e12
end
grid on;
xlim([10, 150]);
legend('r: V_{ex} = 1 mV', 'c: V_{ex} = 1 mV', ...
       'r: V_{ex} = 9 mV', 'c: V_{ex} = 9 mV', ...
       'r: V_{ex} = 20 mV', 'c: V_{ex} = 20 mV');
legend('boxoff');
xlabel('frequency (Hz)');
ylabel('V_{std}/V_{ex}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switch to first measurement configuration (configTG = "top gate only")
configTG.columns  = {'T (K)', 'B (T)', ...
                     'X_{offbal} (V)', 'Y_{offbal} (V)', ...
                     'C/C_{std}', 'C_{loss}/C_{std}', ...
                     'V_{tg} (V)', 'I_{tg} (A)'}; % BG floating
configTG.channels = {@get_probe_temperature, 'AMI430.B', ...
                     'ZI.1X', 'ZI.1Y', ... 
                     'n/a', 'n/a', ...        
                     'K2400.V', 'K2400.I'};
configTG.Vex_amplitude_channel = 'ZI.OUT1AMP1';
configTG.Vex_range_channel     = 'ZI.OUT1RANGE';
configTG.Vstd_amplitude_channel= 'ZI.OUT2AMP2';
configTG.Vstd_range_channel    = 'ZI.OUT2RANGE';
configTG.Vstd_phase_channel    = 'ZI.2phase';
configTG.time_constant_channel = 'ZI.1timeConstant';
configTG.frequency_channel     = 'ZI.1freq';
configTG.Rchip                 = 100e6;
configTG.Vex                   = 9e-3;
configTG.Vstd_range            = 10e-3;
configTG.Xcol = 3;
configTG.Ycol = 4;
configTG.Ccol = 5;
configTG.Lcol = 6;
transistorON;

[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 1e7, 121, configTG, 'log_scale', true);
save('freq_sweep_after_006_configTG_Vex9mV.mat', 'freq', 'Vc0Vex', 'Vr0Vex');

%% test balancing at high freq (prior sweep indicated that Vstd_range was too small for Vex = 9 mV)
smset(configTG.frequency_channel, 10037);
smset(configTG.time_constant_channel, 30e-3);
[balance_matrix, Vc0Vex, Vr0Vex, Cex] = balance_capacitance_bridge(configTG, 'Vex', 8e-3); % 8 mV works

%% try Vtg sweep
configTG.balance_matrix  = balance_matrix;
configTG.Vex             = 8e-3; % new
Vtgcol                   = 7;
Itgcol                   = 8;
configTG.plot_fields     = {configTG.Ccol, configTG.Lcol, Itgcol};
configTG.limit_condition = [Itgcol, 1e-9];

capacitance_bias_sweep(7, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', 0, +1.0, 101, Vtgcol, configTG, 'interval', 1);

%% continue positive
capacitance_bias_sweep(8, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', +1.0, +2.0, 101, Vtgcol, configTG, 'interval', 1);

%% continue positive
capacitance_bias_sweep(9, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', +2.0, +3.0, 101, Vtgcol, configTG, 'interval', 1);

%% no significant leakage yet; let's see negative side before going further
capacitance_bias_sweep(10, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', 0.0, -1.0, 101, Vtgcol, configTG, 'interval', 1);

%% more negative
capacitance_bias_sweep(11, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', -1.0, -2.0, 101, Vtgcol, configTG, 'interval', 1);

%% more negative
capacitance_bias_sweep(12, 'Vtg_sweep_51mK_10.037kHz_Vex8mV', -2.0, -3.0, 101, Vtgcol, configTG, 'interval', 1);
% there doesn't seem to be any change as Vtg is varied from +/-3V
% did the TG lead become disconnected?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switch to BG-only configuration (configBG = "back gate only")
configBG.columns  = {'T (K)', 'B (T)', ...
                     'X_{offbal} (V)', 'Y_{offbal} (V)', ...
                     'C/C_{std}', 'C_{loss}/C_{std}', ...
                     'V_{bg} (V)', 'I_{bg} (A)', 'V_{amp} (V)'}; % TG grounded
configBG.channels = {@get_probe_temperature, 'AMI430.B', ...
                     'ZI.1X', 'ZI.1Y', ... 
                     'n/a', 'n/a', ...        
                     'K2400.V', 'K2400.I', 'K2000.V'};
configBG.Vex_amplitude_channel = 'ZI.OUT1AMP1';
configBG.Vex_range_channel     = 'ZI.OUT1RANGE';
configBG.Vstd_amplitude_channel= 'ZI.OUT2AMP2';
configBG.Vstd_range_channel    = 'ZI.OUT2RANGE';
configBG.Vstd_phase_channel    = 'ZI.2phase';
configBG.time_constant_channel = 'ZI.1timeConstant';
configBG.frequency_channel     = 'ZI.1freq';
configBG.Rchip                 = 100e6;
configBG.Vex                   = 8e-3;
configBG.Vstd_range            = 10e-3;
configBG.Xcol = 3;
configBG.Ycol = 4;
configBG.Ccol = 5;
configBG.Lcol = 6;
transistorON;

[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(10, 1e7, 121, configBG, 'log_scale', true);
save('freq_sweep_after_012_configBG_Vex8mV.mat', 'freq', 'Vc0Vex', 'Vr0Vex');

%% test balancing at high freq (prior sweep indicated that Vstd_range was too small for Vex = 9 mV)
smset(configBG.frequency_channel, 10037);
smset(configBG.time_constant_channel, 30e-3);
[balance_matrix, Vc0Vex, Vr0Vex, Cex] = balance_capacitance_bridge(configBG, 'Vex', 8e-3);

%% try Vbg sweep
configBG.balance_matrix  = balance_matrix;
configBG.Vex             = 8e-3;
Vtgcol                   = 7;
Itgcol                   = 8;
Vampcol                  = 9;
configBG.plot_fields     = {configBG.Ccol, configBG.Lcol, Itgcol, Vampcol};
configBG.limit_condition = [Itgcol, 1e-9];

capacitance_bias_sweep(13, 'Vbg_sweep_51mK_10.037kHz_Vex8mV', 0, +1.0, 101, Vtgcol, configBG, 'interval', 1);
% 1 nA leakage at +0.15 V

%% negative side
capacitance_bias_sweep(14, 'Vbg_sweep_51mK_10.037kHz_Vex8mV', 0.15, -1.0, 101, Vtgcol, configBG, 'interval', 1);
% BG is not going to work either; way too much leakage and a quickly
% changing transistor voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switch to TG/BG-float configuration (configNG = "no gate")
configNG.columns  = {'T (K)', 'B (T)', ...
                     'X_{offbal} (V)', 'Y_{offbal} (V)', ...
                     'C/C_{std}', 'C_{loss}/C_{std}', ...
                     'V_{amp} (V)'}; % TG, BG floating (no DC source connected, only AC excitation on TG)
configNG.channels = {@get_probe_temperature, 'AMI430.B', ...
                     'ZI.1X', 'ZI.1Y', ... 
                     'n/a', 'n/a', ...        
                     'K2000.V'};
configNG.Vex_amplitude_channel = 'ZI.OUT1AMP1';
configNG.Vex_range_channel     = 'ZI.OUT1RANGE';
configNG.Vstd_amplitude_channel= 'ZI.OUT2AMP2';
configNG.Vstd_range_channel    = 'ZI.OUT2RANGE';
configNG.Vstd_phase_channel    = 'ZI.2phase';
configNG.time_constant_channel = 'ZI.1timeConstant';
configNG.frequency_channel     = 'ZI.1freq';
configNG.Rchip                 = 100e6;
configNG.Vex                   = 8e-3;
configNG.Vstd_range            = 10e-3;
configNG.Xcol = 3;
configNG.Ycol = 4;
configNG.Ccol = 5;
configNG.Lcol = 6;
transistorON;

smset(configNG.frequency_channel, 10037);
smset(configNG.time_constant_channel, 300e-3);
[balance_matrix, Vc0Vex, Vr0Vex, Cex] = balance_capacitance_bridge(configNG);
save('balance_matrix_before_015.mat', 'balance_matrix', 'Vc0Vex', 'Vr0Vex', 'Cex');

%% sweep field to see if still measuring capacitance despite gates
Bcol = 2;
smset({'AMI430.rate1', 'AMI430.rate2'}, 0.04); % 10.8 hr back and forth
field_sweep(15, 'field_sweep_configNG_51mK_10.037kHz_Vex8mV', 13, Bcol, configNG, 'plot_fields', {configNG.Xcol, configNG.Ycol});
field_sweep(16, 'field_sweep_configNG_51mK_10.037kHz_Vex8mV', 0, Bcol, configNG, 'plot_fields', {configNG.Xcol, configNG.Ycol});
% error near the beginning of sweep back due to missing log file for
% temperature; modified code to treat this situation in future

%% convert field sweep data to capacitance data
configNG.balance_matrix = balance_matrix;
fnum = 15;
[C, Closs] = offbal2cap_save(fnum, configNG);
figure();
subplot(211); % vert
plot(readcol(fnum, Bcol), C);
subplot(212);
plot(readcol(fnum, Bcol), Closs);
% from the field sweep, I'm not convinced that the TG is still contacting
% the device; there is some change with field but it doesn't look like
% it did last time in the fridge