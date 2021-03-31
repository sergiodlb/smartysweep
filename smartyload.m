function [data, headers, fname] = smartyload(fnum, varargin)
% load data in smartysweep format
%	fnum            <file number; looks for ###*.dat in usual format>
% ---- optional arguments ----
%	data_directory	<path or subdirectory where file is located; may be relative to working directory>
%   header_lines    <number of lines to ignore at the top of the data file; default=1>
%   columns         <array of column numbers to return subset of data in file>
%   fname_format    <format string for data files, if different from default>

% parameters that may change
default_fname_format	= '%03d_*.dat';
default_header_lines	= 1;
default_data_directory  = [];
default_columns         = []; % empty = all columns

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarInt = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'});

% parsed arguments override config fields
addParameter(parser, 'fname_format', default_fname_format, @(s) ischar(s));
addParameter(parser, 'header_lines', default_header_lines, validScalarInt);
addParameter(parser, 'data_directory', default_data_directory);
addParameter(parser, 'columns', default_columns);

parse(parser, varargin{:});
fname_format = parser.Results.fname_format;
header_lines = parser.Results.header_lines;
data_directory = parser.Results.data_directory;
columns = parser.Results.columns;

% load file
fname = sprintf(fname_format, fnum);
if data_directory; fname = fullfile(data_directory, fname); end % build filepath with wildcards
f = dir(fname);
if size(f,1)<1
    if data_directory
        error('No files were found in ''%s'' with format ''%s''.\nCheck <data_directory> or file number.', data_directory, fname);
    else
        error('No files were found in working directory with format ''%s''.\nAre the files located somewhere else?\nAdd ''data_directory'' option or check file number.', fname);
    end
elseif size(f,1)>1
    error('%g files were found with file number %g:\n%s', size(f,1), fnum, sprintf('%s\n', f(:).name));
end
fname = f.name;
if data_directory; fname = fullfile(data_directory, fname); end % build actual filepath
dstr = importdata(fname, '\t', header_lines);
if ~isempty(columns)
    data = dstr.data(:, columns);
else
    data = dstr.data;
end
for col = 1:size(dstr.textdata, 2)-1 % skip timestamp column
    headers{col} = strtrim(dstr.textdata{1,col+1});
end
return