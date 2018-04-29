function chart_recorder(fnum, froot, interval, Npoints, data_fields, plot_fields)
%% start a simple chart recorder
% modified from gate_sweep by Sergio de la Barrera on 2017-01-31
% further modified by Sergio de la Barrera on 2017-03-03 to allow arbitrary data fields
% uses log files to read probe termperature
% - added graceful exit/close plot handling on 2017-03-22
% - modified to handle both get_temperature_from_log used with BlueFors as well as GPIB temps from PPMS or MagLab on 2017-04-24

% parameters that change
log_dir = 'C:\Users\Hunt Lab\Desktop\BlueFors\logs';
log_value = 'T';% log temperature or therm resistance
lw = 2;% plot line width

% % deal with temperature logging method
% if isa(data_fields.channels{1}, 'function_handle')
%     % for log filename
%     log_value = 'T';% log temperature or therm resistance
%     log_date = datestr(clock, 'yy-mm-dd');
%     log_basename = sprintf('CH9 %s %s.log', log_value, log_date);
%     log_path = fullfile(log_dir, log_date, log_basename);
% end

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

% begin recording
verb = 'complete';
start = clock;
tic;
for ii = 1:Npoints
    % build cell array for logging
    dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

    % read smget values
    for col = 1:length(data_fields.columns)
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
                    ax(ll) = plot(data(:, sf), 'LineWidth', lw);
                    hold all;
                end
                xlabel('data points');
                ylabel(data_fields.columns{plot_fields{kk}(1)});
                hold off;
            end
            grid on;
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

    if interval-toc < 0
        disp(sprintf('interval too short by %f s', interval-toc));
    end
    pause(interval-toc);
    tic;
end

% close file
fclose(fid);
fprintf('*** %s\tchart recording %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
end
