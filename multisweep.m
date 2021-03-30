function multisweep(fnum, froot, Vstarts, Vends, Npoints, Vcols, config, varargin)
% performs linear sweep of multiple bias parameters (e.g. DC voltages)
% uses config input structure, a-la the capacitance scripts
% written by Sergio de la Barrera on Jun 23, 2018
%    Vstarts    <starting values for bias parameter; will fast ramp to this point>
%    Vends      <final points for bias parameter>
%    Npoints    <number of points to sweep, including start and end values>
%  **Vcols		<column #s that will be swept by setting bias values>
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and some optionals which can be overridden by varargs:
%                   Vfastrate               (see below)
%                   interval                (see below)
%                   plot_fields             (see below)
%                   plot_xcol               (see below)
%                   limit_condition         (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    dry_run    (flag to simply diplay parameter itinerary without running anything)
%    Vfastrate  <rate at which to fast sweep to starting point; default ~1 V/s>
%    interval   <minimum time in seconds between data points; default = 0>
%    plot_fields    <cell array of columns to live plot; default = {}>
%    plot_xcol      <column # to use as x-axis for plotting; default = Vcol>
%    quiet      <BOOL to block text output to stdout; default = false>
%    limit_condition    <2-element array containing column # and limit value
%                        at which to stop recording and proceed to next point or cancel sweep; default OFF>
%
% 2019-04-24    - moved filename generation to generate_fname()
%               - enabled separate data directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
% deltaVtolerance         = 0.002; %[V]
% default_Vfastrate       = 1; % sweep rate to starting point; approximately ~[V/s]
default_interval        = []; % sweep as quickly as possible unless specified by user
default_plot_fields     = {};
default_plot_xcol       = Vcols(1);
default_quiet           = false; % block all text output (other than errors) if true
default_limit_condition = [];
default_samples         = 1; % for averaging
% Vrange_channel_suffix   = 'range'; % e.g. "range" in "K2400.range"

% validate input parameters
% over_range = false;
% try
%     % check range if instrument supports it
%     V_range = cell2mat(smget([strtok(config.channels{Vcol}, '.'), '.', Vrange_channel_suffix]));
%     if abs(Vstart) > V_range || abs(Vend) > V_range
%         over_range = true;
%     end
% end
% if over_range
%     error('input voltage elements are larger than the available range!');
% end

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validLimitCondition = @(x) validateattributes(x, {'numeric'}, {'numel', 2});
validFunction = @(x) validateattributes(x, {'function_handle'}, {});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run', 'dry-run'})));

% reset defaults based on config entries
% if isfield(config, 'Vfastrate'); default_Vfastrate = config.Vfastrate; end
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end
if isfield(config, 'plot_xcol'); default_plot_xcol = config.plot_xcol; end
if isfield(config, 'limit_condition'); default_limit_condition = config.limit_condition; end

% parsed arguments override config fields
% addParameter(parser, 'Vfastrate', default_Vfastrate, validScalarPos);
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'plot_xcol', default_plot_xcol, validScalarPos); % can override
addParameter(parser, 'limit_condition', default_limit_condition, validLimitCondition); % can override
addParameter(parser, 'samples', default_samples, validScalarNonNeg);
addParameter(parser, 'call_before_measurement', false, validFunction);
addParameter(parser, 'call_after_measurement', false, validFunction);

parse(parser, varargin{:});
% Vfastrate               = parser.Results.Vfastrate;
interval                = parser.Results.interval;
plot_fields             = parser.Results.plot_fields;
plot_xcol               = parser.Results.plot_xcol;
quiet                   = parser.Results.quiet;
limit_condition         = parser.Results.limit_condition;
samples                 = parser.Results.samples;
call_before_measurement = parser.Results.call_before_measurement;
call_after_measurement  = parser.Results.call_after_measurement;

% create list of voltage values
vs = length(Vstarts);
Vlist = zeros(Npoints, vs);
for n = 1:vs
    Vlist(:,n) = linspace(Vstarts(n), Vends(n), Npoints);
end

% execute dry-run
dry_run = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run
    fprintf('DRY-RUN ONLY\n');
    
    % simply plot scan itinerary
    figure();
    plot(Vlist, '-x');
    xlabel('points');
    legend(config.columns{Vcols});
    legend('boxoff');
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

% % fast ramp to start voltage if needed
% V = cell2mat(smget(config.channels{Vcol}));
% if abs(V-Vstart) > deltaVtolerance
%     fprintf('fast ramping %s to %g V\n', config.columns{Vcol}, Vstart);
%     smset(config.channels{Vcol}, Vstart, Vfastrate);
%     while abs(V-Vstart) > deltaVtolerance
%         fprintf('-');
% %         pause(0.2); % doesn't do as intended anyway since above smset holds execution while ramping
%         V = cell2mat(smget(config.channels{Vcol}));
%     end
%     fprintf('> %g V\n', V);
% end

% pre-allocate data array
columns = 1:length(config.columns);
data = zeros(Npoints, length(columns));
data_read = zeros(samples, length(columns));
for n = 1:vs
    data(:, Vcols(n)) = Vlist(:,n); % skip reading sweep parameter; write set values to save precious time
end

% define list of columns to read and fill manual numeric
% normal_columns = columns(columns~=Vcols);
columns(Vcols) = []; % delete columns contained in list, Vcols
normal_columns = columns;
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
start = clock;
tic;
for ii = 1:Npoints
    % go to next voltage (instantaneous)
    smset(config.channels(Vcols), Vlist(ii,:)); % () return a cell array, which allows simultaneous setting
    
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
    if samples > 1
        for sample = 1:samples
            for col = read_columns
                if isa(config.channels{col}, 'function_handle') % could move this logic out of loop
                    data_read(sample, col) = config.channels{col}(); % call user function instead of smget
                else
        %             data_read(sample, col) = cell2mat(smget(config.channels{col}));
                    try                
                        data_read(sample, col) = cell2mat(smget(config.channels{col}));
                    catch gpiberr
                        cprintf('red', 'Warning: error reading channel %s\n', config.channels{col});
                        data_read(sample, col) = nan;
                    end
                end
            end
        end
        data_read(1, :) = nanmean(data_read, 1); % ignore NaN entries
    else
        for col = read_columns
            if isa(config.channels{col}, 'function_handle') % could move this logic out of loop
                data_read(1, col) = config.channels{col}(); % call user function instead of smget
            else
    %             data_read(1, col) = cell2mat(smget(config.channels{col}));
                try                
                    data_read(1, col) = cell2mat(smget(config.channels{col}));
                catch gpiberr
                    cprintf('red', 'Warning: error reading channel %s\n', config.channels{col});
                    data_read(1, col) = nan;
                end
            end
        end
    end
    data(ii, read_columns) = data_read(1, read_columns);
    
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
    
    % test limit condition, if supplied
    if ~isempty(limit_condition) && abs(data(ii, limit_condition(1))) > limit_condition(2)
        if ~quiet
            fprintf('%s is over limit --> %.4g\n', config.columns{limit_condition(1)}, data(ii, limit_condition(1)));
        end
        if ii < Npoints && abs(Vlist(ii+1)) < abs(Vlist(ii)) % if next point is closer to zero
            continue % skip acquisition and go to next point
        else
            break % end execution
        end
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
                    ax(ll) = plot(data(1:ii, plot_xcol), data(1:ii, sf));
                    hold all;
                end
                xlabel(config.columns{plot_xcol});
                ylabel(config.columns{plot_fields{kk}(1)});
                hold off;
            end
            % create swept parameter plot
            subplot(sp_grid{:}, kk+1);
            ax(ll+1) = plot(data(1:ii, plot_xcol), '-k');
            xlabel('data points');
            ylabel(config.columns{plot_xcol});
            sgtitle(fname, 'interpreter', 'none');
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'XData', data(1:ii, plot_xcol), 'YData', data(1:ii, sf));
                    end
                end
                % update swept parameter plot
                set(ax(ll+1), 'YData', data(1:ii, plot_xcol));
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
fprintf('*** %s\tbias sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


