function E_n = volts2E_n_save(fnum, V1col, V2col, config, varargin)
% convert dual gate voltages to electric field and density and save to file

% parameters that may change
time_fixwidth = 24; % must match format strings below
val_fixwidth  = 12; % must match format strings below
format = 'den_%s'; % format of output filename (including original filename)
default_data_directory = [];

% validate required config fields
required_fields = {'Ecol', 'ncol'};
for field = required_fields
    if ~isfield(config, field)
        error('volts2E_n_save requires <%s> in supplied config', char(field));
    end
end
Ecol = config.Ecol;
ncol = config.ncol;

% check for data directory
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
if isfield(config, 'data_directory'); default_data_directory = config.data_directory; end % reset default based on config entry
addParameter(parser, 'data_directory', default_data_directory); % parsed arguments override config fields
parse(parser, varargin{:});
data_directory = parser.Results.data_directory;

% generate data filename
[~, fname_org] = readcol(fnum, 1, varargin{:});
[path, name, ext] = fileparts(fname_org);
nameext = [name, ext];
fname = fullfile(path, sprintf(format, nameext));

% write header
fid = fopen(fname, 'w');
data_header = sprintf('\t%+12s', config.columns{:}); % must match fixed width above
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header); % must match fixed width above
fclose(fid);

% convert gate voltages to E,n
% voltages = readcol(fnum, [V1col,V2col], 'data_directory', data_directory);
voltages = readcol(fnum, [V1col,V2col], varargin{:});
fields = volts2E_n(voltages(:,1), voltages(:,2), config, varargin{:});
E = fields(:,1);
n = fields(:,2);

% read old file line by line and save edited version to new file
fid_org = fopen(fname_org, 'r');
fid = fopen(fname, 'a');
n_idx0 = (time_fixwidth+1) + (ncol-1)*(val_fixwidth+1); % 24 + 1 and 12 + 1 for fixed width plus delimiter '\t' character for each column
E_idx0 = (time_fixwidth+1) + (Ecol-1)*(val_fixwidth+1);
line = fgets(fid_org); % skip header line
ii = 0;
while line ~= -1
    ii = ii + 1;
    line = fgetl(fid_org);
    if ischar(line)
%         fprintf(sprintf('%s\n', line)); % display old data line on stdout
        if n_idx0 > length(line)
            fmt = '\t%12g';
            idx = n_idx0 + (0:12);
        else
            fmt = '%12g';
            idx = n_idx0 + (1:12);
        end
        line(idx) = sprintf(fmt, n(ii));
        if E_idx0 > length(line)
            fmt = '\t%12g';
            idx = E_idx0 + (0:12);
        else
            fmt = '%12g';
            idx = E_idx0 + (1:12);
        end
        line(idx) = sprintf(fmt, E(ii));
        line = sprintf('%s\n', line);
%         fprintf(line); % display new data line on stdout
        fprintf(fid, line);
    end
end
fprintf('saved to --> %s\n', fname);
fclose(fid_org);
fclose(fid);
end

