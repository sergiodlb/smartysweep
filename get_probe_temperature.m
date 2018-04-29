function [ temperature, date, time ] = get_probe_temperature()
% read temperature data from last line of BlueFors logfile

% parameters that change
probe_channel = 9;  % channel of probe thermometer set in LS372
log_value = 'T';    % log temperature or thermometer resistance
log_dir = 'C:\Users\Hunt Lab\Desktop\BlueFors\logs';

% build log path, source for temperature data (has midnight handling; skips current read if log file does not exist)
log_date = datestr(clock, 'yy-mm-dd');
log_basename = sprintf('CH%d %s %s.log', probe_channel, log_value, log_date);
log_path = fullfile(log_dir, log_date, log_basename);

% test existence of log file
if exist(log_path, 'file') ~= 2
    fprintf('temperature log file does not exist yet---is it close to midnight?\nwill attempt to use temperature from just before midnight...\n');
    log_date = datestr(datetime('yesterday'), 'yy-mm-dd'); % use last value from day before if log file does not exist yet
    log_basename = sprintf('CH%d %s %s.log', probe_channel, log_value, log_date);
    log_path = fullfile(log_dir, log_date, log_basename);
    if exist(log_path, 'file') ~= 2
        temperature = nan;
        date = nan;
        time = nan;
        return
    end
end

% access log file
logf = fopen(log_path);
line = fgetl(logf);
while ischar(line)
    tail = line;
    line = fgetl(logf);
end
tail = strsplit(tail, ',');
date = strtrim(tail{1});
time = strtrim(tail{2});
temperature = sscanf(tail{3}, '%f');
fclose(logf);
end

