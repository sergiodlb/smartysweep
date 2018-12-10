%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% example generic configuration for capacitance measurement
config.columns  = {'frequency (Hz)','X_{offbal} (V)',   'Y_{offbal} (V)',   'C/C_{std}',    'C_{loss}/C_{std}'};
config.channels = {'ZI.1freq',      'ZI.1X',            'ZI.1Y',            '',             ''};
config.Vex_amplitude_channel = 'ZI.OUT1AMP1';
config.Vex_range_channel     = 'ZI.OUT1RANGE';
config.Vstd_amplitude_channel= 'ZI.OUT2AMP2';
config.Vstd_range_channel    = 'ZI.OUT2RANGE';
config.Vstd_phase_channel    = 'ZI.2phase';
config.time_constant_channel = 'ZI.1timeConstant';
config.frequency_channel     = 'ZI.1freq';
config.Vex                   = 1e-3;
config.Vstd_range            = 10e-3;
config.fcol = 1;
config.Xcol = 2;
config.Ycol = 3;
config.Ccol = 4;
config.Lcol = 5;

freq_sweep(0, 'test', 10, 100, 5, config.fcol, config, 'plot_fields', {[config.Ccol,config.Lcol]});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% example transistor test
transistor.columns  = {'V_{dd} (V)',	'V_g (V)',	'V_{ds} (V)'};
transistor.channels = {'Vdd',          'Vg',       'Vd'};
transistor.plot_fields = {3};

smset('Vdd', 0.3, 1);
smset('Vg', 0, 1);

Vg_col = 2;
bias_sweep(1, 'transistor_test_Vdd0.3', 0, -0.3, 301, Vg_col, transistor);
transistorOFF;
% Vg* = -0.150 V

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% example measurement of standard capacitor
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

[Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(1, 200, 31, configSTD, 'Vex', 20e-3, 'Vstd_range', 100e-3, 'RC_time', 0.1)
save('measure_Cstd_Vex20mV_35K.mat', 'freq', 'Vc0Vex', 'Vr0Vex');
% Cstd ~ 8.269 pF