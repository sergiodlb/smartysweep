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

freq_sweep(0, 'test', 10, 100, 5, config.fcol, config, 'plot_fields', {2})