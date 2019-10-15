function transistorON(config, varargin)
    % defaults
    Vdd_channel = 'Vdd';
    Vg_channel  = 'Vg';
    Vd_channel  = 'Vd';
    Vdd = 0;
    Vg  = 0;
    
    % deal with optional config
    if nargin > 0
        parser = inputParser;
        parser.KeepUnmatched = true; % other args ignored

        % reset defaults based on config entries
        if isfield(config, 'Vdd_channel'); Vdd_channel = config.Vdd_channel; end
        if isfield(config, 'Vg_channel');  Vg_channel  = config.Vg_channel;  end
        if isfield(config, 'Vd_channel');  Vd_channel  = config.Vd_channel;  end
        if isfield(config, 'Vdd'); Vdd = config.Vdd; end
        if isfield(config, 'Vg');  Vg  = config.Vg;  end

        % parsed arguments override config fields
        addParameter(parser, 'Vdd_channel', Vdd_channel);
        addParameter(parser, 'Vg_channel',  Vg_channel);
        addParameter(parser, 'Vd_channel',  Vd_channel);
        addParameter(parser, 'Vdd', Vdd);
        addParameter(parser, 'Vg', Vg);

        parse(parser, varargin{:});
        Vdd_channel = parser.Results.Vdd_channel;
        Vg_channel  = parser.Results.Vg_channel;
        Vd_channel  = parser.Results.Vd_channel;
        Vdd = parser.Results.Vdd;
        Vg  = parser.Results.Vg;
    end
    
    % set Vg first to limit current through HEMT
    smset(Vg_channel, Vg, 1);
    pause(0.1);
    % then apply drain voltage
    smset(Vdd_channel, Vdd, 1);
    pause(0.1);
    fprintf('transistor ON\namplifier --> %.4g V\n', cell2mat(smget(Vd_channel)));
end