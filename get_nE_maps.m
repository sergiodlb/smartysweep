function [data, n, E, headers, DATA] = get_nE_maps(fnums, config, varargin)
% Converts dual-gate 2D map data to n-E map data and returns 3D array
% Returns:
%   data    [Npoints x Nfiles x Ncolumns]
%   n       [Npoints x Nfiles]
%   E       [Npoints x Nfiles]
%   headers {Ncolumns}

default_data_directory = '';
default_rows = [];
default_cols = [];
default_nshift = false;
fin_format  = '%03d_*.dat';
fout_format = 'nE_maps_%03d-%03d.mat';

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored

% reset defaults based on config entries
if isfield(config, 'data_directory'); default_data_directory = config.data_directory; end

% parsed arguments override config fields
addParameter(parser, 'data_directory', default_data_directory); % can override
addParameter(parser, 'rows', default_rows); % for interpolation
addParameter(parser, 'cols', default_cols); % for interpolation
addParameter(parser, 'nshift', default_nshift); % only affects interpolated data
parse(parser, varargin{:});
data_directory = parser.Results.data_directory;
rows = parser.Results.rows;
cols = parser.Results.cols;
nshift = parser.Results.nshift;

% validate required config fields
required_fields = {'Vbg_col', 'Vtg_col'};
for field = required_fields
    if ~isfield(config, field)
        st = dbstack;
        error('%s(...) requires <%s> in supplied config', st.name, char(field));
    end
end

args = {'data_directory', data_directory, 'fname_format', fin_format};
N = 0;
for fnum = fnums
    fprintf('loading %g...\n', fnum);
    [data, headers, fname] = smartyload(fnum, args{:});

    % convert to E,n
    f = volts2E_n(data(:,config.Vbg_col), data(:,config.Vtg_col), config);
    E = f(:,1); %[V/nm]
    n = f(:,2); %[nm^-2]

    if nshift
        n = n + nshift;
    end
    
    % add to data structure
    N = N + 1;
    DATA(N).data	= data;
    DATA(N).E       = E;
    DATA(N).n       = n;
    DATA(N).headers = headers;
    DATA(N).fnum	= fnum;
    DATA(N).fname	= fname;
end

% gridify the data
try
    % reshape as nd-array, if possible
    data	= cat(3, DATA.data);
    E       = cat(2, DATA.E);
    n       = cat(2, DATA.n);
    data = permute(data, [1,3,2]);
    fprintf('reshaped data\n');
    % If n,E is not on a regular grid, use, for example:
    %   im = pcolor(n, E, data(:,:,4)); set(im, 'edgecolor', 'none');
    % Note that reshaped data may need to be flipped in one or both
    % directions!
catch
    % interpolate irregular-length files
    fprintf('interpolating...');
    data_list	= cat(1, DATA.data);
    E_list      = cat(1, DATA.E);
    n_list      = cat(1, DATA.n);
    
    % assuming E-fixed, n-changing along each sweep
    if isempty(rows)
        rows = max(arrayfun(@(s) size(s.E,1), DATA));
    end
    if isempty(cols)
        cols = length(fnums);
    end
    [n,E] = meshgrid(linspace(min(n_list), max(n_list), rows), ...
                     linspace(min(E_list), max(E_list), cols));
    n = n';
    E = E';

    n_data_columns = size(data_list, 2);
    data = zeros(rows, cols, n_data_columns);
    for N = 1:n_data_columns % loop over columns from original data files
        data(:,:,N) = griddata(n_list, E_list, data_list(:,N), n, E);
    end
    fprintf(' done\n');
    
    % When plotting, remove NaN values with:
    %   im = imagesc(data); set(im, 'alphadata', ~isnan(data));

end

% save collated data
% vars = {'DATA', 'EE', 'NN', 'VBG', 'VTG'};
%     save(fout_name, vars{:});

headers = DATA(1).headers; % return cell array of headers from first file
return

