function flag = cosweep(fnum, froot, V1, V2, V1col, V2col, config, varargin)
% performs a simultaneous bias sweep of two parameters
% based on capactiance_cosweep
% this function written by Sergio de la Barrera on Aug 01, 2018
%    V1         <array of values to sweep for parameter #1>
%    V2         <array of values to sweep for parameter #2; MUST BE EQUAL LENGTH TO V1>
%    V1col		<column # to set V1 values>
%    V2col		<column # to set V2 values>
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and some optionals which can be overridden by varargs:
%                   interval                (see below)
%                   plot_fields             (see below)
%                   limit_condition         (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    dry_run    (flag to simply diplay parameter itinerary without running anything)
%    interval   <minimum time in seconds between data points; default = 0>
%    plot_fields    <cell array of columns to live plot; default = {}>
%    quiet      <BOOL to block text output to stdout; default = false>
%    limit_condition    <2-element array containing column # and limit value
%                        at which to stop recording and proceed to next point or cancel sweep; default OFF>
%
% V1, V2 MUST BE EQUAL LENGTH
% NO AUTO RAMPING
% PARAMETERS ARE SET INSTANTANEOUSLY ON EACH LOOP
% 2018-04-16    - many changes to speed up execution
%               - removed auto wait/settle time based on time constant
%               - moved read_column logic outside of main loop
%               - enabled manual setting of numeric values for columns
%               - added inputParser for optional arguments
%               - added limit_condition option which breaks execution of
%                 main loop on specified data_field rising above limit
% 2018-04-24    - vast reduction in number of positional arguments; now
%                 using fields in "config" structure to convey capacitance
%                 settings
%               - added optional arguments for basic execution options with
%                 ability to override options in config structure by
%                 choosing name-value pair as optional vararg to function
% 2019-04-24    - moved filename generation to generate_fname()
%               - enabled separate data directory
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
validFunction = @(x) validateattributes(x, {'function_handle'}, {});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run', 'dry-run'})));

% reset defaults based on config entries
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end
if isfield(config, 'limit_condition'); default_limit_condition = config.limit_condition; end

% parsed arguments override config fields
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'limit_condition', default_limit_condition, validLimitCondition); % can override
addParameter(parser, 'call_before_measurement', false, validFunction);
addParameter(parser, 'call_after_measurement', false, validFunction);

parse(parser, varargin{:});
interval                = parser.Results.interval;
plot_fields             = parser.Results.plot_fields;
quiet                   = parser.Results.quiet;
limit_condition         = parser.Results.limit_condition;
call_before_measurement = parser.Results.call_before_measurement;
call_after_measurement  = parser.Results.call_after_measurement;

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
fname = generate_fname(fnum, froot, config, varargin{:});

% write header
fid = fopen(fname, 'a');
data_header = sprintf('\t%+12s', config.columns{:});
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header);

% begin sweeping
h = []; % create empty variable for figure handles (for testing whether plot exists later)
verb = 'complete';
flag = 1; % default exit flag
start = clock;
tic;
for ii = 1:length(V1)
    % go to next parameter pair (instantaneous)
    smset({config.channels{V1col}, config.channels{V2col}}, [V1(ii), V2(ii)]);
    
    % call specified function before measurement
    if isa(call_before_measurement, 'function_handle')
        out = call_before_measurement(config, 'return_values', true, varargin{:});
        if isfield(out, 'columns')
            data(ii, out.columns) = out.values;
        end
    end
    
    % build cell array for logging
    dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

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
    
    % call specified function after measurement
    if isa(call_after_measurement, 'function_handle')
        out = call_after_measurement(config, 'return_values', true, varargin{:});
        if isfield(out, 'columns')
            data(ii, out.columns) = out.values;
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
            flag = 0; % issue warning
            break % end execution
        end
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
        if isempty(h)
            % create figure first time
            h = figure();
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
fprintf('*** %s\tco-sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


