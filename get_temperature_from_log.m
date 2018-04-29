function [ temperature, date, time ] = get_temperature_from_log( log_path )
% read temperature data from last line of BlueFors logfile
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

