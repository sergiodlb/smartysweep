function transistorOFF(config, varargin)
    % defaults
    Vdd_channel = 'Vdd';
    Vg_channel  = 'Vg';
    Vd_channel  = 'Vd';
    
    % deal with optional config
    if nargin > 0
        parser = inputParser;
        parser.KeepUnmatched = true; % other args ignored

        % reset defaults based on config entries
        if isfield(config, 'Vdd_channel'); Vdd_channel = config.Vdd_channel; end
        if isfield(config, 'Vg_channel');  Vg_channel  = config.Vg_channel;  end
        if isfield(config, 'Vd_channel');  Vd_channel  = config.Vd_channel;  end

        % parsed arguments override config fields
        addParameter(parser, 'Vdd_channel', Vdd_channel);
        addParameter(parser, 'Vg_channel',  Vg_channel);
        addParameter(parser, 'Vd_channel',  Vd_channel);

        parse(parser, varargin{:});
        Vdd_channel = parser.Results.Vdd_channel;
        Vg_channel  = parser.Results.Vg_channel;
        Vd_channel  = parser.Results.Vd_channel;
    end
    
    smset({Vdd_channel, Vg_channel}, 0, 1);
    fprintf('transistor OFF\namplifier --> %.4g V\n', cell2mat(smget(Vd_channel)));
end