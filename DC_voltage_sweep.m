function DC_voltage_sweep(fnum, froot, Vstart, Vend, interval, Npoints, Vcol, data_fields, plot_fields, xcol, limit_cond)
%% perform a DC bias voltage sweep
% - modified from gate_leakage_test by Sergio de la Barrera on 2017-01-28
% - further modified from gate_sweep by Sergio de la Barrera on 2017-03-01
% - added limit_cond optional argument for short circuit condition on 2017-03-04
% - uses log files to read probe termperature
% - voltage remains on final value
%
%   fnum		<file number>
%   froot		'filename'
%   Vstart		<first voltage of sweep; will fast ramp to this value at start>
%   Vend		<final voltage value; will stop after>
%   interval	<time between data points in seconds>
%   Npoints		<number of data points sweep; interval*Npoints + fastramptime = total time for sweep>
% **Vcol		<column number of voltage for sweeping; contains special measure channel defined in data_fields; PUTTING THE WRONG COLUMN HERE WILL CAUSE SCRIPT TO SMSET ON THE WRONG INSTRUMENT!!!>
%   data_fields <structure with two elements: data_fields.columns = {'header 1', 'header 2', ...} and data_fields.channels = {'INST.CHANNEL1', 'INST.CHANNEL2', ...}>
%   plot_fields {<column>, <numbers>, <to plot>, ..., [<can>, <put>, <in array>, <for same plot>]}
%   xcol        (optional) <column to use as x-axis for plotting>
%   limit_cond  (optional) [<column to track>, <limiting value; will send to zero and stop if above this limit>] 
%
% 2017-03-23    added graceful exit/close plot handling on 
% 2017-04-25    modified to handle both get_temperature_from_log used with
%               BlueFors as well as GPIB temps from PPMS or MagLab
% 2017-09-07    modified to skip read/get of Vcol if read operation is
%               determined not to be supported (does not know state of V)

% parameters that change
Vfastrate = 10;
deltaVtolerance = 0.0001; %[V]
log_dir = 'C:\Users\Hunt Lab\Desktop\BlueFors\logs';
log_value = 'T';% log temperature or therm resistance
% CHECK RANGE OF KEITHLEY BEFORE PROCEEDING!!

% deal with optional plotting argument
if nargin < 10
    xcol = Vcol;
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

% initiate data collection
Vlist = linspace(Vstart, Vend, Npoints);
try
    V = cell2mat(smget(data_fields.channels{Vcol}));
    can_read = true;
    if abs(V-Vstart) > deltaVtolerance
        % fast ramp to start voltage if needed
        disp(sprintf('fast ramping DC bias to %g V', Vstart));
        smset(data_fields.channels{Vcol}, Vstart, Vfastrate);
        while abs(V-Vstart) > deltaVtolerance
            fprintf('-');
%             pause(0.2); % doesn't do as intended anyway since above smset holds execution while ramping
            V = cell2mat(smget(data_fields.channels{Vcol}));
        end
        fprintf('> %g V\n', V);
    end
catch geterr
    % skip if voltage source does not allow read/get
    if strcmp(geterr.message, 'Operation not supported')
        disp(sprintf('Skipping voltage read (not supported for %s)', data_fields.channels{Vcol}));
        can_read = false;
    else
        rethrow(geterr);
    end
end

% pre-allocate data array
columns = 1:length(data_fields.columns);
data = zeros(Npoints, length(columns));
data(:, Vcol) = Vlist; % skip reading sweep parameter; write set values to save precious time

% define list of columns to read and fill manual numeric
normal_columns = columns(columns~=Vcol);
n = 1;
for col = normal_columns
    channel = data_fields.channels{col};
    if isnumeric(channel)
        data(:, col) = channel; % must be scalar or array of appropriate length
        data_fields.columns{col} = ['*', data_fields.columns{col}];
    elseif ~isempty(data_fields.channels{col}) && ~strcmp(data_fields.channels{col}, 'n/a')
        read_columns(n) = col;
        n = n + 1;
    end
end

% generate data filename
fname = sprintf('%03.f_%s.dat', fnum, froot);
while exist(fname, 'file') == 2
%     error('filename exists already');
    fnum = fnum + 1;
    disp(sprintf('*** %s exists already, trying %d', fname, fnum));
    fname = sprintf('%03.f_%s.dat', fnum, froot);
end

% write header
fid = fopen(fname, 'a');
data_header = sprintf('\t%+12s', data_fields.columns{:});
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header);

% begin sweeping
verb = 'complete';
over_limit = false;
start = clock;
tic;
for ii = 1:Npoints
    % go to next voltage (instantaneous)
    smset(data_fields.channels{Vcol}, Vlist(ii));

    % build cell array for logging
    dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

    % read smget values
    for col = read_columns
        if isa(data_fields.channels{col}, 'function_handle')
            % build log path, source for temperature data (has midnight handling; skips current read if log file does not exist)
            log_date = datestr(clock, 'yy-mm-dd');
            log_basename = sprintf('CH9 %s %s.log', log_value, log_date);
            log_path = fullfile(log_dir, log_date, log_basename);
            if exist(log_path, 'file') ~= 2
                disp('temperature log file does not exist yet---is it close to midnight?\twill attempt to use previous temperature...');
            else
                data(ii, col) = data_fields.channels{col}(log_path);
            end
        else
            data(ii, col) = cell2mat(smget(data_fields.channels{col}));
        end
    end
    
    % test limit condition, if supplied
    if (nargin > 10) && (abs(data(ii, limit_cond(1))) > limit_cond(2))
        % set gate voltage to zero
        fprintf('LIMITING VALUE REACHED\nsending %s --> ', data_fields.columns{Vcol});
        smset(data_fields.channels{Vcol}, 0);
        if can_read
            fprintf('%.3g\n', cell2mat(smget(data_fields.channels{Vcol})));
        end
        over_limit = true;
    end

    % record
    data_row = num2cell(data(ii, :));
    data_str = sprintf('\t%12g', data_row{:});
    status = sprintf('%-24s%s\n', dt, data_str);
    fprintf(fid, status);
    fprintf(status);

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
                    ax(ll) = plot(data(1:ii, xcol), data(1:ii, sf));
                    hold all;
                end
                xlabel(data_fields.columns{xcol});
                ylabel(data_fields.columns{plot_fields{kk}(1)});
                hold off;
            end
            % create DC bias plot
            subplot(sp_grid{:}, kk+1);
            ax(ll+1) = plot(data(1:ii, Vcol), '-k');
            xlabel('data points');
            ylabel(data_fields.columns{Vcol});
        else
            % update existing plots with new data
            ll = 0;
            try
                for kk = 1:length(plot_fields)
                    for sf = plot_fields{kk}
                        ll = ll + 1;
                        set(ax(ll), 'XData', data(1:ii, xcol), 'YData', data(1:ii, sf));
                    end
                end
                % update DC bias plot
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

    if over_limit
        error('%s OVER SET LIMIT', data_fields.channels{limit_cond(1)});
    elseif interval-toc < 0
        disp(sprintf('interval too short by %f s', interval-toc));
    end
    pause(interval-toc);
    tic;
end

% close file
fclose(fid);
fprintf('*** %s\tDC voltage bias sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end


