function [C, Closs] = offbal2cap_save(fnum, config, varargin)
% convert off-balance voltages to capacitance and loss and save to file

% parameters that may change
time_fixwidth = 24; % must match format strings below
val_fixwidth  = 12; % must match format strings below

% validate required config fields
required_fields = {'Ccol', 'Lcol'};
for field = required_fields
    if ~isfield(config, field)
        error('offbal2cap_save requires <%s> in supplied config', char(field));
    end
end
Ccol = config.Ccol;
Lcol = config.Lcol;

% generate data filename
[ignore, fname_org] = readcol(fnum, 1);
fname = sprintf('cap_%s', fname_org);

% write header
fid = fopen(fname, 'w');
data_header = sprintf('\t%+12s', config.columns{:}); % must match fixed width above
fprintf(fid, '%-24s%s\n', '#Timestamp', data_header); % must match fixed width above
fclose(fid);

% compute capacitance and loss
[C, Closs] = offbal2cap(fnum, config, varargin{:});

% read old file line by line and save edited version to new file
fid_org = fopen(fname_org, 'r');
fid = fopen(fname, 'a');
C_idx0 = (time_fixwidth+1) + (Ccol-1)*(val_fixwidth+1); % 24 + 1 and 12 + 1 for fixed width plus delimiter '\t' character for each column
Closs_idx0 = (time_fixwidth+1) + (Lcol-1)*(val_fixwidth+1);
line = fgets(fid_org); % skip header line
ii = 0;
while line ~= -1
    ii = ii + 1;
    line = fgetl(fid_org);
    if ischar(line)
%         fprintf(sprintf('%s\n', line)); % display old data line on stdout
        if C_idx0 > length(line)
            fmt = '\t%12g';
            idx = C_idx0 + (0:12);
        else
            fmt = '%12g';
            idx = C_idx0 + (1:12);
        end
        line(idx) = sprintf(fmt, C(ii));
        if Closs_idx0 > length(line)
            fmt = '\t%12g';
            idx = Closs_idx0 + (0:12);
        else
            fmt = '%12g';
            idx = Closs_idx0 + (1:12);
        end
        line(idx) = sprintf(fmt, Closs(ii));
        line = sprintf('%s\n', line);
%         fprintf(line); % display new data line on stdout
        fprintf(fid, line);
    end
end
fprintf('saved to --> %s\n', fname);
fclose(fid_org);
fclose(fid);
end

