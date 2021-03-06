function E_n = volts2E_n(V1, V2, config, varargin)
% for a capacitor with balance point on the middle plate, DC voltages on
% the middle and bottom plates, convert average electric field and electron
% density to V1, V2
e           = 1.602e-19; %[C]
ep0         = 8.854e-21; %[F/nm]
default_ep  = 4;
default_gate_mode = 'bgtg';
default_Vg  = 0;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});

% reset defaults based on config entries
if isfield(config, 'ep'); default_ep = config.ep; end
if isfield(config, 'gate_mode'); default_gate_mode = config.gate_mode; end
if isfield(config, 'Vg'); Vg = config.Vg; else Vg = default_Vg; end

% parsed arguments override config fields
addParameter(parser, 'ep', default_ep, validScalarPos); % can override
addParameter(parser, 'gate_mode', default_gate_mode, @(x) any(validatestring(x, {'bgdc','tgdc','bgtg'}))); % can override
parse(parser, varargin{:});
gate_mode = parser.Results.gate_mode;

% validate required config fields (some required by balance_capacitance_bridge)
required_fields = {'dbg', 'dtg'};
% if strcmp(gate_mode, 'capacitance')
%     required_fields{end+1} = 'Vg';
% end
for field = required_fields
    if ~isfield(config, field)
        error('E_n_2volts requires <%s> in %s mode', char(field), gate_mode);
    end
end

% define local variables
ep	= parser.Results.ep;
dbg = config.dbg;
dtg = config.dtg;

% transform to E, n
if strcmp(gate_mode, 'bgdc') % V1 = Vbg, V2 = V_{middle-layer}
%     Vg	= config.Vg;
    E = 0.5*((V1 - V2)/dbg - (Vg - V2)/dtg);
    n = ep*ep0*((V1 - V2)/dbg + (Vg - V2)/dtg)/e;
elseif strcmp(gate_mode, 'tgdc') % V1 = Vtg, V2 = V_{middle-layer}
    E = 0.5*((Vg - V2)/dbg - (V1 - V2)/dtg);
    n = ep*ep0*((Vg - V2)/dbg + (V1 - V2)/dtg)/e;
else % V1 = bg, V2 = Vtg
    E = 0.5*((V1-Vg)/dbg - (V2-Vg)/dtg);
    n = ep*ep0*((V1-Vg)/dbg + (V2-Vg)/dtg)/e;
end
E_n = [E, n];
end