function freq_sweep(fnum, froot, fstart, fend, Npoints, fcol, config, varargin)
%% perform a frequency sweep
% written by Sergio de la Barrera on Mar 3, 2017
%    fnum		<file number>
%    froot		'filename'
%    fstart		<first frequency value; will jump to this value at start>
%    fend		<final frequency value; will stop after>
%    Npoints    <number of data points sweep; total time for sweep also affected by changing time constant and settle time>
%  **fcol		<column # of frequency channel for sweeping; used by smset>
%    config     structure containing:
%                   channels = {'INST.CHANNEL1','INST.CHANNEL2', ...}
%                   columns  = {'header 1',     'header 2',      ...}
%               and some optionals which can be overridden by varargs:
%                   interval                (see below)
%                   plot_fields             (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    log_scale  <BOOL toggle log scale for frequency sweep; default = true>
%    interval   <minimum time in seconds between data points; default = 0>
%    tc_mult    <prefactor to set time constant based on freq; tc = tc_mult/freq; default = 3>
%    plot_fields    <cell array of columns to live plot; default = {}>
%                   e.g. {<column>, <numbers>, <to plot>, ..., [<can>, <put>, <in array>, <for same plot>]}
%    quiet      <BOOL to block text output to stdout; default = false>
%
% 2018-07-06    - updated to use config structure and key-val optional
%                 parameters
% 2018-07-20    - added call_before_measurement and call_after_measurement
%                 optional parameters that will execute a specified
%                 function call and store any returned values in the data
%                 columns specified in config + called function
% 2019-04-24    - moved filename generation to generate_fname()
%               - enabled separate data directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
default_log_scale       = true; % toggle log x-scale and log spacing of freq points
default_interval        = []; % sweep as quickly as possible unless specified by user
default_set_tc          = true; % choose whether or not to change tc based on frequency
default_tc_mult         = 3; % prefactor to set time constant based on freq; tc = tc_mult/freq
default_wait_for_lock   = false; % if given an instrument channel, will query and wait for confirmation of reference lock from lock-in amplifier
default_plot_fields     = {};
default_quiet           = false; % block all text output (other than errors) if true

% validate required config fields (some required by balance_capacitance_bridge)
required_fields = {'time_constant_channel'};
for field = required_fields
    if ~isfield(config, field)
        error('freq_sweep requires <%s> in supplied config', char(field));
    end
end

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validFunction = @(x) validateattributes(x, {'function_handle'}, {});

addParameter(parser, 'log_scale', default_log_scale);
addParameter(parser, 'set_tc', default_set_tc);
addParameter(parser, 'tc_mult', default_tc_mult);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'call_before_measurement', false, validFunction);
addParameter(parser, 'call_after_measurement', false, validFunction);

% reset defaults based on config entries
if isfield(config, 'wait_for_lock'); default_wait_for_lock = config.wait_for_lock; end
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end

% parsed arguments override config fields
addParameter(parser, 'wait_for_lock', default_wait_for_lock, @ischar);
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override

parse(parser, varargin{:});
log_scale               = parser.Results.log_scale;
interval                = parser.Results.interval;
set_tc                  = parser.Results.set_tc;
tc_mult                 = parser.Results.tc_mult;
wait_for_lock           = parser.Results.wait_for_lock;
plot_fields             = parser.Results.plot_fields;
quiet                   = parser.Results.quiet;
call_before_measurement = parser.Results.call_before_measurement;
call_after_measurement  = parser.Results.call_after_measurement;

% choose subplot layout
np = length(plot_fields);
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

% list of data columns
columns = 1:length(config.columns);

% generate data filename
fname = generate_fname(fnum, froot, config, varargin{:});

% write header
fid = fopen(fname, 'a');
data_header = sprintf('\t%+12s', config.columns{:});
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header);

% initiate data collection
freq = logspace(log10(fstart), log10(fend), Npoints);

% toggle frequency log scale
if log_scale
    freq = logspace(log10(fstart), log10(fend), Npoints);
else
    freq = linspace(fstart, fend, Npoints);
end

% pre-allocate measurement arrays
data = zeros(Npoints, length(columns));
data(:, fcol) = freq; % skip reading sweep parameter; write set values to save precious time

% define list of columns to read and fill manual numeric
normal_columns = columns(columns~=fcol);
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

% begin sweeping
verb = 'complete';
start = clock;
tic;
for ii = 1:Npoints
    % go to next frequency (instantaneous)
    smset(config.channels{fcol}, freq(ii));
    if set_tc; smset(config.time_constant_channel, tc_mult/freq(ii)); end % set time_constant based on freq
    pause(tc_mult/freq(ii)); % wait to settle
    if wait_for_lock % if given a channel for lock-in reference lock (could be any boolean query)
        while ~cell2mat(smget(wait_for_lock)) % if unlocked, wait
            fprintf('waiting for reference lock\n');
            pause(tc_mult/freq(ii)); % wait to settle
        end
    end
    
    % call specified function before measurement
    if isa(call_before_measurement, 'function_handle')
        out = call_before_measurement(config, 'return_values', true, varargin{:});
        if isfield(out, 'columns')
            data(ii, out.columns) = out.values;
        end
    end

    % build cell array for logging
    dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

    % read smget values
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
                    ax(ll) = plot(data(1:ii, fcol), data(1:ii, sf));
                    hold all;
                end
                xlabel(config.columns{fcol});
                ylabel(config.columns{plot_fields{kk}(1)});
                hold off;
                if log_scale
                    set(gca, 'XScale', 'log');
                end
                sgtitle(fname, 'interpreter', 'none');
            end
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'XData', data(1:ii, fcol), 'YData', data(1:ii, sf));
                    end
                end
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
fprintf('*** %s\tfrequency sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


