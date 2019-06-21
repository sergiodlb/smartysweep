function fname = generate_fname(fnum, froot, config, varargin)
% simple function used by all measurement routines to check data directory,
% file numbers, and generate unique filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deal with optionals
fstring                 = '%03.f_%s.dat';
default_data_directory  = [];

parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
if isfield(config, 'data_directory'); default_data_directory = config.data_directory; end % reset defaults based on config entries
addParameter(parser, 'data_directory', default_data_directory); % parsed arguments override config fields
parse(parser, varargin{:});
data_directory = parser.Results.data_directory;

% create data directory if necessary
if data_directory
    data_directory = fullfile(data_directory); % normalize the file path
    if exist(data_directory, 'dir') ~= 7 % code for existing directory
        mkdir(data_directory);
        fprintf('created data directory --> ''%s''\n', data_directory);
    end
end

% create unique filename (based on entire name, not just number)
fname = sprintf(fstring, fnum, froot);
if data_directory; fname = fullfile(data_directory, fname); end
while exist(fname, 'file') == 2
    fnum = fnum + 1;
    fprintf('*** %s exists already, trying %d\n', fname, fnum);
    fname = sprintf(fstring, fnum, froot);
    if data_directory; fname = fullfile(data_directory, fname); end
end

return

