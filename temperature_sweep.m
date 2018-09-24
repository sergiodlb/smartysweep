function temperature_sweep(fnum, froot, Tset, config, varargin)
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
%
% 2018-06-23 re-written to use config input structure and key-value
%            optional argument pairs
% 2018-06-24 added sweep_direction logic to end recording as soon as T is
%            beyond set point, regardless of sweep direction
% 2018-07-20    - added call_before_measurement and call_after_measurement
%                 optional parameters that will execute a specified
%                 function call and store any returned values in the data
%                 columns specified in config + called function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOES NOT SET TEMPERATURE; MANUALLY SET SYSTEM TO CONDENSE, WARM, OR COOL

% parameters that change
deltaTtolerance          = 0.001; % absolute temp difference in Kelvin
default_interval         = []; % measure as frequently as possible unless specified by user
default_maxruntime       = inf; % only quit on reaching Tset unless a max timeout condition is given [hours]
default_plot_fields      = {};
default_quiet            = false; % block all text output (other than errors) if true
default_play_sound       = false;
default_send_email_after = {}; % do not send an email unless address(es) given by user

% validate required config fields
required_fields = {'Tcol'};
for field = required_fields
    if ~isfield(config, field)
        error('temperature_sweep requires <%s> in supplied config', char(field));
    end
end
Tcol = config.Tcol;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
validFunction = @(x) validateattributes(x, {'function_handle'}, {});

% reset defaults based on config entries
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end

% parsed arguments override config fields
addParameter(parser, 'interval', default_interval, validScalarNonNeg); % can override
addParameter(parser, 'maxruntime', default_maxruntime, validScalarPos);
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'play_sound', default_play_sound);
addParameter(parser, 'send_email_after', default_send_email_after, @iscell);
addParameter(parser, 'call_before_measurement', false, validFunction);
addParameter(parser, 'call_after_measurement', false, validFunction);

parse(parser, varargin{:});
interval                = parser.Results.interval;
maxruntime              = parser.Results.maxruntime;
plot_fields             = parser.Results.plot_fields;
quiet                   = parser.Results.quiet;
play_sound              = parser.Results.play_sound;
send_email_after        = parser.Results.send_email_after;
call_before_measurement = parser.Results.call_before_measurement;
call_after_measurement  = parser.Results.call_after_measurement;

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
columns = 1:length(config.columns);

% generate data filename
fname = sprintf('%03.f_%s.dat', fnum, froot);
while exist(fname, 'file') == 2
    fnum = fnum + 1;
    disp(sprintf('*** %s exists already, trying %d', fname, fnum));
    fname = sprintf('%03.f_%s.dat', fnum, froot);
end

% write header
fid = fopen(fname, 'a');
data_header = sprintf('\t%+12s', config.columns{:});
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
if isa(config.channels{Tcol}, 'function_handle')
    T = config.channels{Tcol}();
else
    T = cell2mat(smget(config.channels{Tcol}));
end
sweep_direction = sign(Tset-T);

% begin data collection
verb = 'complete';
ii = 0;
last_Ttime = '';
status = '';
start = clock;
tic;
try
    while (abs(T-Tset) > deltaTtolerance && (T-Tset)*sweep_direction < 0) && etime(clock, start) < maxruntime*3600 && ~isnan(T)
        if isa(config.channels{Tcol}, 'function_handle')
            [T, Tdate, Ttime] = config.channels{Tcol}();
            if ~strcmp(Ttime, last_Ttime)% && T ~= 0
                % since the lastest entry is indeed new, record the new data point
                last_Ttime = Ttime;
                record = true;
            end
        else
%             T = cell2mat(smget(config.channels{Tcol}));
            try                
                T = cell2mat(smget(config.channels{Tcol}));
                record = true;
            catch gpiberr
                cprintf('red', 'Warning 1: error reading channel %s\n', config.channels{Tcol});
                T = nan;
                record = false;
            end
%             record = true;
        end
        if record
            ii = ii + 1;
            
            % call specified function before measurement
            if isa(call_before_measurement, 'function_handle')
                out = call_before_measurement(config, 'return_values', true, varargin{:});
                if isfield(out, 'columns')
                    data(ii, out.columns) = out.values;
                end
                % get temperature again in case significant time has passed
                if isa(config.channels{Tcol}, 'function_handle')
                    [T, Tdate, Ttime] = config.channels{Tcol}();
                    if ~strcmp(Ttime, last_Ttime)
                        last_Ttime = Ttime;
                    end
                else
%                     T = cell2mat(smget(config.channels{Tcol}));
                    try                
                        T = cell2mat(smget(config.channels{Tcol}));
                    catch gpiberr
                        cprintf('red', 'Warning 2: error reading channel %s\n', config.channels{Tcol});
                        T = nan;
                    end
                end
            end
            
            % build cell array for logging
            dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');
            if ~isnan(T) && ~isempty(T) % still giving errors?
                data(ii, Tcol) = T;
            else % assign nan if empty due to GPIB error
                data(ii, Tcol) = nan;
            end
            
            % read other smget channels
            for col = columns(columns~=Tcol)
                channel = config.channels{col};
                if isa(channel, 'function_handle')
                    data(ii, col) = channel(); % call user function instead of smget
                elseif isnumeric(channel) % TO-DO: this is outside of loop in other functions
                    data(ii, col) = channel; % must be scalar
                    if ii == 1
                        config.columns{col} = ['*', config.columns{col}];
                    end
                elseif isempty(channel) || strcmp(channel, 'n/a') % do nothing, defaults entries to zero (possible conflict with call_before_measurement otherwise)
%                     data(ii, col) = 0; % skip columns with empty channel name or 'n/a'
                else
%                     data(ii, col) = cell2mat(smget(config.channels{col}));
                    try                
                        data(ii, col) = cell2mat(smget(config.channels{col}));
                    catch gpiberr
                        cprintf('red', 'Warning 3: error reading channel %s\n', config.channels{col});
                        data(ii, col) = nan;
                    end
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
                        % plot selected columns vs temperature
                        subplot(sp_grid{:}, kk);
                        for sf = plot_fields{kk}
                            ll = ll + 1;
                            ax(ll) = plot(data(:, Tcol), data(:, sf));
                            hold all;
                        end
                        xlabel(config.columns{Tcol});
                        ylabel(config.columns{plot_fields{kk}(1)});
                        hold off;
                    end
                    % create temperature plot
                    subplot(sp_grid{:}, kk+1);
                    ax(ll+1) = plot(data(:, Tcol), '-k');
                    xlabel('data points');
                    ylabel(config.columns{Tcol});
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
        
        if interval
            if ~quiet && interval-toc < 0
                disp(sprintf('interval too short by %f s', interval-toc));
            end
            pause(interval-toc);
            tic;
        end
    end

    % reset
    message = sprintf('*** %s\ttemperature sweep %s\n%s', datestr(clock, 'mmm dd HH:MMPM'), verb, status);
    subject = sprintf('data collection %s, T = %g', verb, T);
    disp(message);
catch err
%     try
%         message = sprintf('*** error has occurred\n%s\n%s\n%s', err.identifier, err.message, status);
%     catch err2
%         message = sprintf('*** error has occurred\n%s\n%s\n%s', err.identifier, err.message);
%     end
    message = sprintf('*** error has occurred\n%s\n%s\n%s', err.identifier, err.message);
    err.stack
    T
    subject = sprintf('measurement error has occurred');
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
