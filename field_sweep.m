function field_sweep(fnum, froot, Bset, Bcol, config, varargin)
%% measure during a field sweep
% modified from gate_sweep and RT_measurement by Sergio de la Barrera on 2017-01-30
% uses log files to read probe termperature
% waits for field change to start and end
% 2017-03-24 added smset functionality for initiating ramp to field target
%            correct ramp rates must be set by user beforehand
%            *** note that above relies on Bcol to be accurate!!!
% 2017-04-24 modified to handle both get_temperature_from_log used with BlueFors as well as GPIB temps from PPMS or MagLab 
% 2018-04-05 modified to skip 'n/a' channels (yields zero)
% 2018-04-20 implemented numeric channel input and get_probe_temperature()
% 2018-04-25 major changes including:
%            - reduction in number of positional arguments; now
%              using fields in "config" structure to convey many settings
%            - added optional arguments for basic execution options with
%              ability to override options in config structure by
%              choosing name-value pair as optional vararg to function
% 2018-06-24 added sweep_direction logic to end recording as soon as T is
%            beyond set point, regardless of sweep direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
deltaBtolerance          = 0.0005; %[B]
default_interval         = []; % measure as frequently as possible unless specified by user
default_plot_fields      = {};
default_quiet            = false; % block all text output (other than errors) if true
default_set_field        = true;
default_play_sound       = false;
default_notify_when      = 0.95; % play audible sound when this fraction complete
default_send_email_after = {}; % do not send an email unless address(es) given by user

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarFrac = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', '<=', 1.0});

% reset defaults based on config entries
if isfield(config, 'interval'); default_interval = config.interval; end
if isfield(config, 'plot_fields'); default_plot_fields = config.plot_fields; end

% parsed arguments override config fields
addParameter(parser, 'interval', default_interval, validScalarNonNeg); % can override
addParameter(parser, 'plot_fields', default_plot_fields, @iscell); % can override
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'set_field', default_set_field);
addParameter(parser, 'play_sound', default_play_sound);
addParameter(parser, 'notify_when', default_notify_when, validScalarFrac);
addParameter(parser, 'send_email_after', default_send_email_after, @iscell);

parse(parser, varargin{:});
interval                = parser.Results.interval;
plot_fields             = parser.Results.plot_fields;
quiet                   = parser.Results.quiet;
set_field               = parser.Results.set_field;
play_sound              = parser.Results.play_sound;
notify_when             = parser.Results.notify_when;
send_email_after        = parser.Results.send_email_after;

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
    fprintf('*** %s exists already, trying %d\n', fname, fnum);
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

% get initial field
readB = cell2mat(smget(config.channels{Bcol}));
B = readB;
Binit = readB;
sweep_direction = sign(Bset-Binit);

if set_field
    % set field target
    smset(config.channels{Bcol}, Bset);
end

% begin measuring
verb = 'complete';
oblivious = true;
ii = 0;
start = clock;
tic;
while abs(readB-Bset) > deltaBtolerance && (readB-Bset)*sweep_direction < 0
    % read field value
    readB = cell2mat(smget(config.channels{Bcol}));

    % wait for field change before recording
    if abs(readB-B(end)) > deltaBtolerance
        % field value is new, continue
        ii = ii + 1;
        B(ii) = readB;
        
        % audible progress indicator
        if play_sound && oblivious && (readB-Binit)/(Bset-Binit) > notify_when
            oblivious = false;
            load splat;
            soundsc(y, Fs);
        end

        % build cell array for logging
        dt = datestr(clock, 'yyyy-mm-ddTHH:MM:SS.FFF');

        % read smget values
        for col = columns
            channel = config.channels{col};
            if isa(channel, 'function_handle')
                data(ii, col) = channel(); % call user function instead of smget
            elseif isnumeric(channel)
                data(ii, col) = channel; % must be scalar
                if ii == 1
                    config.columns{col} = ['*', config.columns{col}];
                end               
            elseif isempty(channel) || strcmp(channel, 'n/a')
                data(ii, col) = 0; % skip columns with empty channel name or 'n/a'
            else
                data(ii, col) = cell2mat(smget(config.channels{col}));                    
            end
        end

        % record
        data_row = num2cell(data(ii, :));
        data_str = sprintf('\t%12g', data_row{:});
        status = sprintf('%-24s%s\n', dt, data_str);
        fprintf(fid, status);
        if ~quiet
            fprintf(status);
        end;

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
                        ax(ll) = plot(data(:, Bcol), data(:, sf));
                        hold all;
                    end
                    xlabel(config.columns{Bcol});
                    ylabel(config.columns{plot_fields{kk}(1)});
                    hold off;
                end
                % create field plot
                subplot(sp_grid{:}, kk+1);
                ax(ll+1) = plot(data(:, Bcol), '-k');
                xlabel('data points');
                ylabel(config.columns{Bcol});
            else
                % update existing plots with new data
                ll = 0;
                try
                    for kk = 1:length(plot_fields)
                        for sf = plot_fields{kk}
                            ll = ll + 1;
                            set(ax(ll), 'XData', data(:, Bcol), 'YData', data(:, sf));
                        end
                    end
                    % update temperature plot
                    set(ax(ll+1), 'YData', data(:, Bcol));
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
end

% close file
fclose(fid);
subject = sprintf('BlueFors field sweep %s, B = %g T', verb, readB);
message = sprintf('%s\tfield sweep %s\n', datestr(clock, 'mmm dd HH:MMPM'), verb);
fprintf('*** %s', message);

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
