function out = measure_capacitance(config, varargin)
% converts measured off-balance voltages to capacitance
%   Xcol            <column of offbal voltage X-component in file>
%   Ycol            <column of offbal voltage Y-component in file>
%   balance_matrix  <output of balance_capacitance_bridge(...)>
%   Vex             (optional)<excitation voltage; req for real capacitance>
%   Cstd            (optional)<standard capacitance; req for real capacitance>
%   return_values   (optional)<toggles measurement mode (default = false), which will return column numbers and Xoffbal, Yoffbal, C, Closs values if called by another measurement function>
%                   e.g.
%                     if return_values
%                         out.columns = [config.Xcol, config.Ycol, config.Ccol, config.Lcol]; % tell parent function where to put data values
%                         out.values  = [Xoffbal, Yoffbal, Cex, Closs]; % return data values to be stored in the specified columns
%                     else
%                         out = [Cex, Closs]; % simply return measured capacitance and loss
%                     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that may change
default_Vex          = []; % throw error if empty
default_Cstd         = 1; % if standard capacitance is unknown

% validate required config fields
required_fields = {'balance_matrix', 'Xcol', 'Ycol'};
for field = required_fields
    if ~isfield(config, field)
        error('calculate_capacitance requires <%s> in supplied config', char(field));
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
addParameter(parser, 'return_values', false); % true if called by a measurement script (e.g. to compute capacitance from live off-balance measurement)

parse(parser, varargin{:});
Cstd = parser.Results.Cstd;

% validate excitation voltage
Vex = parser.Results.Vex;
if isempty(Vex)
    error('calculate_capacitance requires Vex in supplied config or as optional argument');
end

% check for additional arguments needed if returning data values
return_values = parser.Results.return_values;
if return_values
    required_fields = {'Ccol', 'Lcol'};
    for field = required_fields
        if ~isfield(config, field)
            error('calculate_capacitance requires <%s> in supplied config in order to return data for measurement', char(field));
        end
    end
end 

% unpack balance matrix
balance_matrix = num2cell(balance_matrix);
[Kc1, Kc2, Kr1, Kr2, Vc0, Vr0] = balance_matrix{:};

% compute all the necessaries
Xoffbal = cell2mat(smget(config.channels{Xcol})); % new value will replace any prev measured one in live measurement
Yoffbal = cell2mat(smget(config.channels{Ycol})); % new value will replace any prev measured one in live measurement
L1prime = Xoffbal;
L2prime = Yoffbal;
Vr0prime = Vr0 + (Kc2 * L1prime - Kc1 * L2prime) / (Kc1 * Kr2 - Kr1 * Kc2);
Vc0prime = Vc0 + (Kr1 * L2prime - Kr2 * L1prime) / (Kc1 * Kr2 - Kr1 * Kc2);

% calculate capacitance
Cex = Cstd * Vc0prime / Vex; % edit on 1/8/2019 for tuning antisymetric capacitance (can be negative)
Closs = Cstd * Vr0prime / Vex;

% handle data output for live measurement versus direct call by user
if return_values
    out.columns = [config.Xcol, config.Ycol, config.Ccol, config.Lcol]; % tell parent function where to put data values
    out.values  = [Xoffbal, Yoffbal, Cex, Closs]; % return data values to be stored in the specified columns
else
    out = [Cex, Closs]; % simply return measured capacitance and loss
end
return