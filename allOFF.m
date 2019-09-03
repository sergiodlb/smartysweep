function allOFF(config)

% defaults
Vdd_channel  = 'Vdd';
Vg_channel   = 'Vg';
Vex_channel  = [];
Vstd_channel = [];
Vinv_channel = [];
Vbg_channel  = [];
Vtg_channel  = [];
Vdc_channel  = [];

if nargin > 0
    if isfield(config, 'Vdd_channel'); Vdd_channel = config.Vdd_channel; end
    if isfield(config, 'Vg_channel');  Vg_channel  = config.Vg_channel;  end
    if isfield(config, 'Vbg_col');     Vbg_channel = config.channels{config.Vbg_col}; end
    if isfield(config, 'Vtg_col');     Vtg_channel = config.channels{config.Vtg_col}; end
    if isfield(config, 'Vdc_col');     Vdc_channel = config.channels{config.Vdc_col}; end
    if isfield(config, 'Vex_amplitude_channel');  Vex_channel  = config.Vex_amplitude_channel; end
    if isfield(config, 'Vstd_amplitude_channel'); Vstd_channel = config.Vstd_amplitude_channel; end
    if isfield(config, 'Vinv_channel'); Vinv_channel = config.Vinv_channel; end
end
list = {Vdd_channel, Vg_channel, Vbg_channel, Vtg_channel, Vdc_channel, Vex_channel, Vstd_channel, Vinv_channel};
list = list(~cellfun('isempty',list)); % remove empty elements

fprintf('all OFF..!\n');
try
    smset(list, 0, 0.5);
    fprintf('%s\n', sprintf('%12s\t', list{:}));
    fprintf('%s\n', sprintf('%12.g\t', cell2mat(smget(list))));
catch
    smset(list, 0);
    fprintf('%s set to zero\n', sprintf('%12s, ', list{:}));
end
end