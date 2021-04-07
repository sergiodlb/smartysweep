function Vbg_Vdc = E_n2volts(E, n, config, varargin)
% for a capacitor with balance point on the middle plate, DC voltages on
% the middle and bottom plates, convert average electric field and electron
% density to Vbg, Vdc
e           = 1.602e-19; %[C]
ep0         = 8.854e-21; %[F/nm]
default_ep  = 4;
default_gate_mode = 'bgtg';
default_Vg  = 0;

% deal with optional arguments
parser = inputParser;
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
Cb  = ep*ep0/dbg;
Ct  = ep*ep0/dtg;

% transform to Vbg, Vdc
if strcmp(gate_mode, 'bgdc')
%     Vg	= config.Vg;
    Vbg = +0.5*((2*ep*ep0*E + e*n)/Cb + (2*ep*ep0*E - e*n)/Ct) + Vg;
    Vdc = +0.5*((2*ep*ep0*E - e*n)/Ct) + Vg;
elseif strcmp(gate_mode, 'tgdc') % Vbg means Vtg here
    Vbg = -0.5*((2*ep*ep0*E + e*n)/Cb + (2*ep*ep0*E - e*n)/Ct) + Vg;
    Vdc = -0.5*((2*ep*ep0*E + e*n)/Cb) + Vg;
else  % Vdc means Vtg here
    Vbg = +0.5*((2*ep*ep0*E + e*n)/Cb) + Vg;
    Vdc = -0.5*((2*ep*ep0*E - e*n)/Ct) + Vg; % Vdc --> Vtg
end
Vbg_Vdc = [Vbg, Vdc]; % Vdc --> Vtg if 'bgtg' mode
end
