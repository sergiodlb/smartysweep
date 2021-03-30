function chart_recorder(fnum, froot, config, varargin)
%% start a simple chart recorder
% uses config input structure, a-la the capacitance scripts
% written by Sergio de la Barrera on Jan 01, 2017
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and some optionals which can be overridden by varargs:
%                   interval                (see below)
%                   plot_fields             (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    interval   <minimum time in seconds between data points; default = 0>
%    Npoints    <finite number of data points to record; default = inf>
%    plot_fields    <cell array of columns to live plot; default = {}>
%    quiet      <BOOL to block text output to stdout; default = false>
%
% ---- change log
% 2017-01-31 first version, modified from gate_sweep script
% 2017-03-03 modified to allow arbitrary data fields
% 2017-03-22 uses log files to read probe termperature
% 			 added graceful exit/close plot handling
% 2017-04-24 modified to handle both get_temperature_from_log used with
%			 BlueFors as well as GPIB temps from PPMS or MagLab
% 2018-06-26 heavily modified to match config input structure and key-value
%			 optional parameters of capacitance scripts
% 2018-07-20    - added call_before_measurement and call_after_measurement
%                 optional parameters that will execute a specified
%                 function call and store any returned values in the data
%                 columns specified in config + called function
% 2019-06-21    - moved filename generation to generate_fname()
%               - enabled separate data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
default_interval        = []; % sweep as quickly as possible unless specified by user
default_Npoints         = inf; % record forever, unless a finite number of points is given
default_plot_fields     = {};
default_quiet           = false; % block all text output (other than errors) if true

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validFunction = @(x) validateattributes(x, {'function_handle'}, {});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run', 'dry-run'})));

% reset defaults based on config entries
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end

% parsed arguments override config fields
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'Npoints', default_Npoints, validScalarPos); % could be improved with integer test
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'call_before_measurement', false, validFunction);
addParameter(parser, 'call_after_measurement', false, validFunction);

parse(parser, varargin{:});
interval                = parser.Results.interval;
Npoints                 = parser.Results.Npoints;
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

% define list of columns to read and fill manual numeric
columns = 1:length(config.columns);
n = 1;
for col = columns
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

% begin recording
verb = 'complete';
start = clock;
tic;
ii = 0;
while ii < Npoints
	ii = ii + 1;

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
                    ax(ll) = plot(data(:, sf));
                    hold all;
                end
                xlabel('data points');
                ylabel(config.columns{plot_fields{kk}(1)});
                hold off;
            end
            sgtitle(fname, 'interpreter', 'none');
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'YData', data(:, sf));
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
            disp(sprintf('interval too short by %f s', interval-toc));
        end
        pause(interval-toc);
        tic;
    end
end

% close file
fclose(fid);
fprintf('*** %s\tchart recording %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end
