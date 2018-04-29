function capacitance_bias_sweep(fnum, froot, Vstart, Vend, Npoints, Vcol, config, varargin)
% performs capacitance measurement while sweeping a bias parameter using 
% prior capacitance bridge balance parameters
% inspired by offbalbiassweep written by Nick Gomon
% this function written by Sergio de la Barrera on Aug 29, 2017
%    Vstart     <starting value for bias parameter; will fast ramp to this point>
%    Vend       <final point for bias parameter>
%    Npoints    <number of points to sweep, including start and end values>
%    Vcol		<column # that will be swept by setting bias values>
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
%                   Vfastrate               (see below)
%                   interval                (see below)
%                   plot_fields             (see below)
%                   calculate_capacitance   (see below)
%                   rebalance               (see below)
%                   Cstd                    (see below)
%                   limit_condition         (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    dry_run    (flag to simply diplay parameter itinerary without running anything)
%    Vfastrate  <rate at which to fast sweep to starting point; default ~1 V/s>
%    interval   <minimum time in seconds between data points; default = 0>
%    plot_fields    <cell array of columns to live plot; default = {}>
%    quiet      <BOOL to block text output to stdout; default = false>
%    calculate_capacitance  <BOOL to compute capacitance from off-balance voltages; default = true>
%    rebalance  <BOOL to rebalance at each sweep point; default = false>
%    Cstd       <standard capacitance; default = 1>
%    limit_condition    <2-element array containing column # and limit value
%                        at which to stop recording and proceed to next point or cancel sweep; default OFF>
%
% 2018-03-26    - revised to use voltage units (a la balance_capacitance_bridge)
%               - added rebalance option (hardcode)
% 2018-03-27    - further corrections to units; essentially replaced
%                 vr0, vr0prime and such with real voltage quantities
%               - capacatitance is now Cstd * abs(Vc0prime / Vex)
% 2018-04-05    - shortened settle time to 3*time_const
% 2018-04-06    - eliminated Gcol and references to conductance, now
%                 requires Lcol and outputs capacitance loss (no freq req)
% 2018-04-16    - many changes to speed up execution
%               - removed auto wait/settle time based on time constant
%               - moved read_column logic outside of main loop
%               - made calculating capacitance (or leaving zero) an option
%               - enabled manual setting of numeric values for columns
% 2018-04-25    - vast reduction in number of positional arguments; now
%                 using fields in "config" structure to convey capacitance
%                 settings
%               - added optional arguments for basic execution options with
%                 ability to override options in config structure by
%                 choosing name-value pair as optional vararg to function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
deltaVtolerance         = 0.0001; %[V]
default_Vfastrate       = 1; % sweep rate to starting point; approximately ~[V/s]
default_interval        = []; % sweep as quickly as possible unless specified by user
default_plot_fields     = {};
default_quiet           = false; % block all text output (other than errors) if true
default_calculate_capacitance = true;
default_rebalance       = false;
default_Cstd            = 1;
default_limit_condition = [];
Vrange_channel_suffix   = 'range'; % e.g. "range" in "K2400.range"

% validate input parameters
try
    % check range if instrument supports it
    V_range = smget([strtok(config.channels{Vcol}, '.'), '.', Vrange_channel_suffix]);
    if abs(Vstart) > V_range || abs(Vend) > V_range
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
if isfield(config, 'Vfastrate'); default_Vfastrate = config.Vfastrate; end
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end
if isfield(config, 'calculate_capacitance'); default_calculate_capacitance = config.calculate_capacitance; end
if isfield(config, 'rebalance'); default_rebalance = config.rebalance; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end
if isfield(config, 'limit_condition'); default_limit_condition = config.limit_condition; end

% parsed arguments override config fields
addParameter(parser, 'Vfastrate', default_Vfastrate, validScalarPos);
addParameter(parser, 'interval', default_interval, validScalarNonNeg);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'calculate_capacitance', default_calculate_capacitance); % can override
addParameter(parser, 'rebalance', default_rebalance); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override
addParameter(parser, 'limit_condition', default_limit_condition, validLimitCondition); % can override

parse(parser, varargin{:});
Vfastrate               = parser.Results.Vfastrate;
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

% create list of voltage values
Vlist = linspace(Vstart, Vend, Npoints);

% execute dry-run
dry_run = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run
    fprintf('DRY-RUN ONLY\n');
    
    % simply plot scan itinerary
    figure();
    plot(Vlist, '-x');
    xlabel('points');
    ylabel(config.columns{Vcol});
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

% fast ramp to start voltage if needed
V = cell2mat(smget(config.channels{Vcol}));
if abs(V-Vstart) > deltaVtolerance
    fprintf('fast ramping %s to %g V\n', config.columns{Vcol}, Vstart);
    smset(config.channels{Vcol}, Vstart, Vfastrate);
    while abs(V-Vstart) > deltaVtolerance
        fprintf('-');
%         pause(0.2); % doesn't do as intended anyway since above smset holds execution while ramping
        V = cell2mat(smget(config.channels{Vcol}));
    end
    fprintf('> %g V\n', V);
end

% pre-allocate data array
columns = 1:length(config.columns);
data = zeros(Npoints, length(columns));
data(:, Vcol) = Vlist; % skip reading sweep parameter; write set values to save precious time

% define list of columns to read and fill manual numeric
normal_columns = columns(columns~=Vcol);
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
for ii = 1:Npoints
    % go to next voltage (instantaneous)
    smset(config.channels{Vcol}, Vlist(ii));
    
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
                    ax(ll) = plot(data(1:ii, Vcol), data(1:ii, sf));
                    hold all;
                end
                xlabel(config.columns{Vcol});
                ylabel(config.columns{plot_fields{kk}(1)});
                hold off;
            end
            % create swept parameter plot
            subplot(sp_grid{:}, kk+1);
            ax(ll+1) = plot(data(1:ii, Vcol), '-k');
            xlabel('data points');
            ylabel(config.columns{Vcol});
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'XData', data(1:ii, Vcol), 'YData', data(1:ii, sf));
                    end
                end
                % update swept parameter plot
                set(ax(ll+1), 'YData', data(1:ii, Vcol));
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
fprintf('*** %s\tcapacitance bias sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


