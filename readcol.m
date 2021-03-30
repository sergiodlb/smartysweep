function [data, fname] = readcol(fnum, cols, varargin)
% short function to read a specific column from data in smartysweep format
%	fnum            <file number; looks for ###*.dat in usual format>
%   cols            <array of column numbers of data in file to return>
% ---- optional arguments ----
%	data_directory	<path or subdirectory where file is located; may be relative to working directory>
%   header_lines    <number of lines to ignore at the top of the data file; default=1>
%   fname_format    <format string for data files, if different from default>

% parameters that may change
default_fname_format	= '%03d_*.dat';
default_header_lines	= 1;
default_data_directory  = [];

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarInt = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'});

% parsed arguments override config fields
addParameter(parser, 'fname_format', default_fname_format, @(s) ischar(s));
addParameter(parser, 'header_lines', default_header_lines, validScalarInt);
addParameter(parser, 'data_directory', default_data_directory);

parse(parser, varargin{:});
fname_format = parser.Results.fname_format;
header_lines = parser.Results.header_lines;
data_directory = parser.Results.data_directory;

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
data = dstr.data(:, cols);
return