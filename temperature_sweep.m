function temperature_sweep(fnum, froot, Tset, interval, maxruntime, data_fields, plot_fields, send_email_after)
%% measure in BlueFors during warmup or cooldown; requires manual change in temperature
% written by Sergio de la Barrera on 2017-01-16
% modified on 2017-01-26
% heavily modified on 2017-02-25
% - converted to function
% - interval is in seconds
% - maxruntime is in hours
% - data fields to record are provided as structure with
%   data_fields.columns containing the names in the header file (and plot labels)
%   and data_fields.channels containing a cell array of the smget channels
% - email is provided as optional argument at end
% - added midnight handling for get_temperature_from_logs; same as others
% - added graceful exit/close plot handling on 2017-03-22
% - modified to handle both get_temperature_from_log used with BlueFors as well as GPIB temps from PPMS or MagLab on 2017-04-24
% 2018-04-05 modified to skip 'n/a' channels (yields zero)

% parameters that change
Tcol = 1;
deltaTtolerance = 0.001; % absolute temp difference in Kelvin
use_resistance = false;
log_dir = 'C:\Users\Hunt Lab\Desktop\BlueFors\logs';
play_sound = true;
% MANUALLY SET SYSTEM TO CONDENSE!
% WRITE DOWN FREQ & TERMINAL CONFIGURATION!

% deal with optional email argument
if nargin < 8
    send_email_after = {};
end

% deal with temperature logging method
if isa(data_fields.channels{1}, 'function_handle')
    % for log filename
    if use_resistance
        log_value = 'R';
    else
        log_value = 'T';
    end
    log_date = datestr(clock, 'yy-mm-dd');
    log_basename = sprintf('CH9 %s %s.log', log_value, log_date);
    log_path = fullfile(log_dir, log_date, log_basename);
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

% list of data columns
columns = 1:length(data_fields.columns);

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

% initialize email notifications
if ~isempty(send_email_after)
    faddress = 'huntlabcmu@gmail.com';

    % obtain and encrypt password
    import javax.crypto.Cipher;

    cipher = Cipher.getInstance('RSA');
    keygen = java.security.KeyPairGenerator.getInstance('RSA');
    keyPair = keygen.genKeyPair();
    cipher.init(Cipher.ENCRYPT_MODE, keyPair.getPrivate());

    plaintextUnicodeVals = uint16(passcode);
    plaintextBytes = typecast(plaintextUnicodeVals, 'int8');
    ciphertext = cipher.doFinal(plaintextBytes)';
end

% get initial temperature
if isa(data_fields.channels{Tcol}, 'function_handle')
    T = data_fields.channels{Tcol}(log_path);
else
    T = cell2mat(smget(data_fields.channels{Tcol}));
end

% begin data collection
verb = 'complete';
ii = 0;
last_Ttime = '';
status = '';
start = clock;
tic;
try
    while abs(T-Tset) > deltaTtolerance && etime(clock, start) < maxruntime*3600
        if isa(data_fields.channels{Tcol}, 'function_handle')
            % build log path, source for temperature data
            log_date = datestr(clock, 'yy-mm-dd');
            log_basename = sprintf('CH9 %s %s.log', log_value, log_date);
            log_path = fullfile(log_dir, log_date, log_basename);
            if exist(log_path, 'file') ~= 2
                disp('temperature log file does not exist yet---is it close to midnight?\twill attempt to use previous temperature...');
            else
                [T, Tdate, Ttime] = get_temperature_from_log(log_path);
            end
            if ~strcmp(Ttime, last_Ttime)% && T ~= 0
                % deal with T > 100 K (TOVER) in BlueFors with LakeShore 372
                if T == 0 && ~use_resistance
                    T = 100;
                end

                % since the lastest entry is indeed new, record the new data point
                last_Ttime = Ttime;
                record = true;
            end
        else
            T = cell2mat(smget(data_fields.channels{Tcol}));
            record = true;
        end
        if record
            ii = ii + 1;
            
            % build cell array for logging
            dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');
            data(ii, Tcol) = T;
            
            % read other smget channels
            for col = columns(columns~=Tcol)
                if strcmp(data_fields.channels{col}, 'n/a') % skip columns with channel name 'n/a'
                    data(ii, col) = 0;
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
                        % plot selected columns vs temperature
                        subplot(sp_grid{:}, kk);
                        for sf = plot_fields{kk}
                            ll = ll + 1;
                            ax(ll) = plot(data(:, Tcol), data(:, sf));
                            hold all;
                        end
                        xlabel(data_fields.columns{Tcol});
                        ylabel(data_fields.columns{plot_fields{kk}(1)});
                        hold off;
                    end
                    % create temperature plot
                    subplot(sp_grid{:}, kk+1);
                    ax(ll+1) = plot(data(:, Tcol), '-k');
                    xlabel('data points');
                    ylabel(data_fields.columns{Tcol});
                else
                    % update existing plots with new data
                    ll = 0;
                    try
                        for kk = 1:length(plot_fields)
                            for sf = plot_fields{kk}
                                ll = ll + 1;
                                set(ax(ll), 'XData', data(:, Tcol), 'YData', data(:, sf));
                            end
                        end
                        % update temperature plot
                        set(ax(ll+1), 'YData', data(:, Tcol));
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
        end
        if interval-toc < 0
            disp(sprintf('interval too short by %f s', interval-toc));
        end
        pause(interval-toc);
        tic;
    end

    % reset
    message = sprintf('*** %s\tR-T measurement %s\n%s', datestr(clock, 'mmm dd HH:MMPM'), verb, status);
    subject = sprintf('BlueFors data collection %s, T = %g', verb, T);
    disp(message);
catch err
    try
        message = sprintf('*** error has occurred\n%s\n%s\n%s', err.identifier, err.message, status);
    catch err2
        message = sprintf('*** error has occurred\n%s\n%s\n%s', err.identifier, err.message);
    end
    subject = sprintf('BlueFors measurement error has occurred');
    disp(message);
end

% close file
fclose(fid);

% audio notification
if play_sound
    load gong;
    soundsc(y, Fs);
end

% send follow-up email
if ~isempty(send_email_after)
    setpref('Internet', 'E_mail', faddress);
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'SMTP_Username', faddress);
    cipher.init(Cipher.DECRYPT_MODE, keyPair.getPublic());
    setpref('Internet', 'SMTP_Password', char(typecast(cipher.doFinal(ciphertext), 'uint16'))');
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth', 'true');
    props.setProperty('mail.smtp.socketFactory.class', ...
                      'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port', '465');
    sendmail(send_email_after, subject, message);
end
