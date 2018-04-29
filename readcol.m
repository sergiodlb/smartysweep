function [data, fname] = readcol(fnum, col)
% short function to read a specific column from data in our typical format
%   fnum            <file number; looks for ###*.dat in usual format>
%   col             <column number of data in file>

% parameters that may change
fname_format = sprintf('%03d*.dat', fnum);
header_lines = 1;

% load file
f = dir(fname_format);
fname = f.name;
dstr = importdata(fname, '\t', header_lines);
data = dstr.data(:, col);
return