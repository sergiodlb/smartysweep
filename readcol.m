function [data, fname] = readcol(fnum, cols, varargin)
% short function to read a specific column from data in our typical format
%   fnum            <file number; looks for ###*.dat in usual format>
%   cols            <array of column numbers of data in file>

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
fname = f.name;
if data_directory; fname = fullfile(data_directory, fname); end % build actual filepath
dstr = importdata(fname, '\t', header_lines);
data = dstr.data(:, cols);
return