function capacitance_cosweep(fnum, froot, V1, V2, V1col, V2col, config, varargin)
% performs capacitance measurement while co-sweeping two parameters
% based on modified capactiance_DC_bias_sweep
% this function written by Sergio de la Barrera on Apr 10, 2018
%    V1         <array of values to sweep for parameter #1>
%    V2         <array of values to sweep for parameter #2; MUST BE EQUAL LENGTH TO V1>
%    V1col		<column # to set V1 values>
%    V2col		<column # to set V2 values>
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and channels specific to capacitance (req to calculate C):
%                   Vex_amplitude_channel   <sm channel of excitation>
%                   Vex_range_channel       <sm channel of Vex range>
%                   Xcol       <column # of off-balance X measurement>
%                   Ycol       <column # of off-balance Y measurement>
%                   Ccol	   <column # to store capacitance>
%                   Lcol	   <column # to store capacitance loss>
%               if not rebalancing, also requires:
%                   balance_matrix  <matrix of cap-bridge balanced parameters: [Kc1, Kc2, Kr1, Kr2, vc0, vr0]>
%               and some optionals which can be overridden by varargs:
%                   interval                (see below)
%                   plot_fields             (see below)
%                   calculate_capacitance   (see below)
%                   rebalance               (see below)
%                   Cstd                    (see below)
%                   limit_condition         (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    dry_run    (flag to simply diplay parameter itinerary without running anything)
%    interval   <minimum time in seconds between data points; default = 0>
%    plot_fields    <cell array of columns to live plot; default = {}>
%    quiet      <BOOL to block text output to stdout; default = false>
%    calculate_capacitance  <BOOL to compute capacitance from off-balance voltages; default = true>
%    rebalance  <BOOL to rebalance at each sweep point; default = false>
%    Cstd       <standard capacitance; default = 1>
%    limit_condition    <2-element array containing column # and limit value
%                        at which to stop recording and proceed to next point or cancel sweep; default OFF>
%
% V1, V2 MUST BE EQUAL LENGTH
% NO AUTO RAMPING
% PARAMETERS ARE SET INSTANTANEOUSLY ON EACH LOOP
% 2018-04-16    - many changes to speed up execution
%               - removed auto wait/settle time based on time constant
%               - moved read_column logic outside of main loop
%               - made calculating capacitance (or leaving zero) an option
%               - enabled manual setting of numeric values for columns
%               - added inputParser for optional arguments; altered Cstd to
%                 use this new paradigm as a name-value pair
%               - added limit_condition option which breaks execution of
%                 main loop on specified data_field rising above limit
% 2018-04-24    - vast reduction in number of positional arguments; now
%                 using fields in "config" structure to convey capacitance
%                 settings
%               - added optional arguments for basic execution options with
%                 ability to override options in config structure by
%                 choosing name-value pair as optional vararg to function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
default_interval        = []; % sweep as quickly as possible unless specified by user
default_plot_fields     = {};
default_quiet           = false; % block all text output (other than errors) if true
default_calculate_capacitance = true;
default_rebalance       = false;
default_Cstd            = 1;
default_limit_condition = [];
Vrange_channel_suffix   = 'range'; % e.g. "range" in "K2400.range"

% validate input parameters
if length(V1) ~= length(V2)
    error('length of V1 and V2 parameters lists must match!');
end
try
    % check range if instrument supports it
    V1_range = smget([strtok(config.channels{V1col}, '.'), '.', Vrange_channel_suffix]);
    V2_range = smget([strtok(config.channels{V2col}, '.'), '.', Vrange_channel_suffix]);
    if any(abs(V1) > V1_range) || any(abs(V2) > V2_range)
        error('input voltage elements are larger than the available range!');
    end
end

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validLimitCondition = @(x) validateattributes(x, {'numeric'}, {'numel', 2});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run'})));

% reset defaults based on config entries
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end
if isfield(config, 'calculate_capacitance'); default_calculate_capacitance = config.calculate_capacitance; end
if isfield(config, 'rebalance'); default_rebalance = config.rebalance; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end
if isfield(config, 'limit_condition'); default_limit_condition = config.limit_condition; end

% parsed arguments override config fields
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'calculate_capacitance', default_calculate_capacitance); % can override
addParameter(parser, 'rebalance', default_rebalance); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override
addParameter(parser, 'limit_condition', default_limit_condition, validLimitCondition); % can override

parse(parser, varargin{:});
interval                = parser.Results.interval;
plot_fields             = parser.Results.plot_fields;
quiet                   = parser.Results.quiet;
calculate_capacitance   = parser.Results.calculate_capacitance;
rebalance               = parser.Results.rebalance;
Cstd                    = parser.Results.Cstd;
limit_condition         = parser.Results.limit_condition;

% validate capacitance_settings
if calculate_capacitance
    required_fields = {'Vex_amplitude_channel', 'Vex_range_channel', 'Xcol', 'Ycol', 'Ccol', 'Lcol'};
    for field = required_fields
        if ~isfield(config, field)
            error('capacitance calculation requires <%s> in supplied config', char(field));
        end
    end
    Xcol = config.Xcol;
    Ycol = config.Ycol;
    Ccol = config.Ccol;
    Lcol = config.Lcol;
    
    if ~rebalance
        if ~isfield(config, 'balance_matrix')
            error('capacitance calculation requires <balance_matrix> in supplied config');
        end
        % unpack balance matrix
        balance_matrix = num2cell(config.balance_matrix);
        [Kc1, Kc2, Kr1, Kr2, Vc0, Vr0] = balance_matrix{:}; % uses actual voltages
    end
    
    % get AC excitation of device under test
    vd = cell2mat(smget(config.Vex_amplitude_channel)); % not true voltage (hence lowercase) --> dimensionless fraction of output voltage scale (ZI only)
    Vex_range = cell2mat(smget(config.Vex_range_channel));
    Vex = vd*Vex_range; % excitation in proper voltage units
end

% execute dry-run
dry_run = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run
    fprintf('DRY-RUN ONLY\n');
    
    % simply plot scan itinerary
    figure();
    plot(V1, V2, '-x');
    xlabel(config.columns{V1col});
    ylabel(config.columns{V2col});
    return % quit execution
end

% choose subplot layout
np = length(plot_fields) + 1;
switch np
    case {1, 2, 3}
        sp_grid = {1, np};
    case {4, 6}
        sp_grid = {2, np/2};
    case 5
        sp_grid = {2, 3};
    case {7, 8, 9}
        sp_grid = {3, 3};
    otherwise
        sp_grid = {1, 1};
end

% pre-allocate data array
columns = 1:length(config.columns);
data = zeros(length(V1), length(columns));
data(:, V1col) = V1; % skip reading sweep parameters; write set values to save precious time
data(:, V2col) = V2;

% define list of columns to read and fill manual numeric
normal_columns = columns(columns~=V1col & columns~=V2col);
n = 1;
for col = normal_columns
    channel = config.channels{col};
    if isnumeric(channel)
        data(:, col) = channel; % must be scalar or array of appropriate length
        config.columns{col} = ['*', config.columns{col}];
    elseif ~isempty(config.channels{col}) && ~strcmp(config.channels{col}, 'n/a')
        read_columns(n) = col;
        n = n + 1;
    end
end

% generate data filename
fname = sprintf('%03.f_%s.dat', fnum, froot);
while exist(fname, 'file') == 2
    fnum = fnum + 1;
    fprintf('*** %s exists already, trying %d\n', fname, fnum);
    fname = sprintf('%03.f_%s.dat', fnum, froot);
end

% write header
fid = fopen(fname, 'a');
data_header = sprintf('\t%+12s', config.columns{:});
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header);

% begin sweeping
verb = 'complete';
start = clock;
tic;
for ii = 1:length(V1)
    % go to next parameter pair (instantaneous)
    smset(config.channels{V1col}, V1(ii));
    smset(config.channels{V2col}, V2(ii));

    % build cell array for logging
    dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

    % read smget values (pre-allocated cols defined before loop)
    for col = read_columns
        if isa(config.channels{col}, 'function_handle') % could move this logic out of loop
            data(ii, col) = config.channels{col}(); % call user function instead of smget
        else
            data(ii, col) = cell2mat(smget(config.channels{col}));
        end
    end
    
    % test limit condition, if supplied
    if ~isempty(limit_condition) && abs(data(ii, limit_condition(1))) > limit_condition(2)
        if ~quiet
            fprintf('%s is over limit --> %.4g\n', config.columns{limit_condition(1)}, data(ii, limit_condition(1)));
        end
        if ii < length(V1) && norm([V1(ii+1), V2(ii+1)]) < norm([V1(ii), V2(ii)]) % if next point is closer to zero
            continue % skip acquisition and go to next point
        else
            break % end execution
        end
    end

    if calculate_capacitance
        if rebalance
            % rebalance at each Vdc point
            [balance_matrix, Vc0Vex, Vr0Vex, Cex] = balance_capacitance_bridge(config, varargin{:});
            Vr0prime = Vr0Vex*Vex;
            Vc0prime = Vc0Vex*Vex;
        else
            % measure off-balance voltage components
            L1prime = 1.005 * sqrt(2) * data(ii, Xcol); % compute L from single measurement directly
            L2prime = 1.005 * sqrt(2) * data(ii, Ycol); % real voltage units
            Vr0prime = Vr0 + (Kc2 * L1prime - Kc1 * L2prime) / (Kc1 * Kr2 - Kr1 * Kc2);
            Vc0prime = Vc0 + (Kr1 * L2prime - Kr2 * L1prime) / (Kc1 * Kr2 - Kr1 * Kc2);
        end

        % calculate device capacitance
        data(ii, Ccol) = Cstd * abs(Vc0prime / Vex);
        data(ii, Lcol) = Cstd * abs(Vr0prime / Vex);
    end

    % record
    data_row = num2cell(data(ii, :));
    data_str = sprintf('\t%12g', data_row{:});
    status = sprintf('%-24s%s\n', dt, data_str);
    fprintf(fid, status);
    if ~quiet
        fprintf(status);
    end

    % plot
    if ~isempty(plot_fields)
        if ii == 1
            % create figure first time
            figure();
            ll = 0;
            for kk = 1:length(plot_fields)
                % plot selected columns vs field
                subplot(sp_grid{:}, kk);
                for sf = plot_fields{kk}
                    ll = ll + 1;
                    ax(ll) = plot(data(1:ii, V1col), data(1:ii, sf));
                    hold all;
                end
                xlabel(config.columns{V1col});
                ylabel(config.columns{plot_fields{kk}(1)});
                hold off;
            end
            % create parameter plot
            subplot(sp_grid{:}, kk+1);
            ax(ll+1) = plot(data(1:ii, V1col), data(1:ii, V2col), '-k');
            xlabel(config.columns{V1col});
            ylabel(config.columns{V2col});
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'XData', data(1:ii, V1col), 'YData', data(1:ii, sf));
                    end
                end
                % update DC bias plot
                set(ax(ll+1), 'XData', data(1:ii, V1col), 'YData', data(1:ii, V2col));
                drawnow;
            catch ploterr
                % exit gracefully if plot closed
                if strcmp(ploterr.identifier, 'MATLAB:class:InvalidHandle')
                    verb = 'canceled';
                    break;
                else
                    rethrow(ploterr);
                end
            end
        end
    end

    if interval
        if ~quiet && interval-toc < 0
            disp(sprintf('interval too short by %f s', toc-interval));
        end
        pause(interval-toc);
        tic;
    end
end

% close file
fclose(fid);
fprintf('*** %s\tcapacitance co-sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


