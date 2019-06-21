function [C, L] = offbal2cap_subtract_spurious(fnum, config, varargin)
% subtract spurious signals (wth scaling) from off-balance voltages and
% converts corrected off-balance voltages to capacitance
%   fnum            <file number; looks for ###*.dat in usual format>
%   Xcol            <column of offbal voltage X-component in file>
%   Ycol            <column of offbal voltage Y-component in file>
%   Xsp_col         <column of spurious signal X-component in file>
%   Ysp_col         <column of spurious signal Y-component in file>
%   balance_matrix  <output of balance_capacitance_bridge(...)>
%   Vex             (optional)<excitation voltage; req for real capacitance>
%   Cstd            (optional)<standard capacitance; req for real capacitance>

% parameters that may change
fname_format         = sprintf('%03d*.dat', fnum);
default_data_directory = [];
default_Vex          = []; % throw error if empty
default_Cstd         = 1; % if standard capacitance is unknown
default_header_lines = 1;
default_Xdata        = [];
default_Ydata        = [];
default_Xsp_col      = [];
default_Ysp_col      = [];
default_xmult        = [];
default_ymult        = [];
default_vppfactor    = 1; % 1.005*sqrt(2) for case where Vex is peak-to-peak voltage; 1 for RMS voltages
default_plot_subtraction = false;
default_save         = false;

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
validScalar = @(x) validateattributes(x, {'numeric'}, {'scalar'});
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validScalarInt = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'});

% reset defaults based on config entries
if isfield(config, 'Vex'); default_Vex = config.Vex; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end
if isfield(config, 'Xsp_col'); default_Xsp_col = config.Xsp_col; end
if isfield(config, 'Ysp_col'); default_Ysp_col = config.Ysp_col; end
if isfield(config, 'data_directory'); default_data_directory = config.data_directory; end % reset default based on config entry

% parsed arguments override config fields
addParameter(parser, 'Vex', default_Vex, validScalarNonNeg); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override
addParameter(parser, 'header_lines', default_header_lines, validScalarInt);
addParameter(parser, 'Xdata', default_Xdata, @(x) isvector(x));
addParameter(parser, 'Ydata', default_Ydata, @(x) isvector(x));
addParameter(parser, 'Xsp_col', default_Xsp_col, validScalarInt);
addParameter(parser, 'Ysp_col', default_Ysp_col, validScalarInt);
addParameter(parser, 'xmult', default_xmult, validScalar);
addParameter(parser, 'ymult', default_ymult, validScalar);
addParameter(parser, 'vppfactor', default_vppfactor, validScalarNonNeg);
addParameter(parser, 'plot_subtraction', default_plot_subtraction);
addParameter(parser, 'save', default_save);
addParameter(parser, 'data_directory', default_data_directory); % parsed arguments override config fields

parse(parser, varargin{:});
Vex = parser.Results.Vex;
Cstd = parser.Results.Cstd;
header_lines = parser.Results.header_lines;
Xdata = parser.Results.Xdata;
Ydata = parser.Results.Ydata;
Xsp_col = parser.Results.Xsp_col;
Ysp_col = parser.Results.Ysp_col;
xmult = parser.Results.xmult;
ymult = parser.Results.ymult;
vppfactor = parser.Results.vppfactor;
plot_subtraction = parser.Results.plot_subtraction;
data_directory = parser.Results.data_directory;
if parser.Results.save
    offbal2cap_fn = @offbal2cap_save;
else
    offbal2cap_fn = @offbal2cap;
end

% validate required parameter fields
required_fields = {'Vex', 'Xs_col', 'Ys_col', 'xmult', 'ymult'};
for field = required_fields
    if isempty(field)
        error('offbal2cap requires <%s> in supplied config or as optional argument', char(field));
    end
end

% load the off-balance data
fprintf('loading %g...\n', fnum);
data = readcol(fnum, [Xcol, Ycol, Xsp_col, Ysp_col], 'data_directory', data_directory);

X = data(:,1);
Y = data(:,2);
Xsp = data(:,3);
Ysp = data(:,4);
x = X - xmult*Xsp;
y = Y - ymult*Ysp;

[C, L] = offbal2cap_fn(fnum, config, 'Xdata', x, 'Ydata', y);
if plot_subtraction
    figure(2);
    subplot(3,2,1); plot(X);
    subplot(3,2,2); plot(Y);
    subplot(3,2,3); plot(Xsp);
    subplot(3,2,4); plot(Ysp);
    subplot(3,2,5); plot(x);
    subplot(3,2,6); plot(y);
end

return