function [Cex, Closs] = offbal2cap(fnum, config, varargin)
% converts off-balance voltages to capacitance
%   fnum            <file number; looks for ###*.dat in usual format>
%   Xcol            <column of offbal voltage X-component in file>
%   Ycol            <column of offbal voltage Y-component in file>
%   balance_matrix  <output of balance_capacitance_bridge(...)>
%   Vex             (optional)<excitation voltage; req for real capacitance>
%   Cstd            (optional)<standard capacitance; req for real capacitance>

% parameters that may change
fname_format         = sprintf('%03d*.dat', fnum);
default_Vex          = []; % throw error if empty
default_Cstd         = 1; % if standard capacitance is unknown
default_header_lines = 1;

% validate required config fields
required_fields = {'balance_matrix', 'Xcol', 'Ycol'};
for field = required_fields
    if ~isfield(config, field)
        error('offbal2cap requires <%s> in supplied config', char(field));
    end
end
balance_matrix = config.balance_matrix;
Xcol = config.Xcol;
Ycol = config.Ycol;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validScalarInt = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'});

% reset defaults based on config entries
if isfield(config, 'Vex'); default_Vex = config.Vex; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end

% parsed arguments override config fields
addParameter(parser, 'Vex', default_Vex, validScalarNonNeg); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override
addParameter(parser, 'header_lines', default_header_lines, validScalarInt);

parse(parser, varargin{:});
Cstd = parser.Results.Cstd;
header_lines = parser.Results.header_lines;

% validate excitation voltage
Vex = parser.Results.Vex;
if isempty(Vex)
    error('offbal2cap requires Vex in supplied config or as optional argument');
end

% load file
f = dir(fname_format);
dstr = importdata(f.name, '\t', header_lines);
data = dstr.data;

% unpack balance matrix
balance_matrix = num2cell(balance_matrix);
[Kc1, Kc2, Kr1, Kr2, Vc0, Vr0] = balance_matrix{:};

% compute all the necessaries
L1prime = 1.005 * sqrt(2) * data(:, Xcol);
L2prime = 1.005 * sqrt(2) * data(:, Ycol);
Vr0prime = Vr0 + (Kc2 * L1prime - Kc1 * L2prime) / (Kc1 * Kr2 - Kr1 * Kc2);
Vc0prime = Vc0 + (Kr1 * L2prime - Kr2 * L1prime) / (Kc1 * Kr2 - Kr1 * Kc2);

% calculate capacitance
Cex = Cstd * abs(Vc0prime / Vex);
Closs = Cstd * abs(Vr0prime / Vex);
return