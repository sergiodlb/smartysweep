function determine_slope(fnums, data_col, config, varargin)
% Plots dual-gate 2D map data to aid in determining ratio of capacitances
default_data_directory = '';
default_cmap = 'hot';
default_marker = 10;
default_one_panel = false;
default_scatter_only = false;
default_draggable_line = true;
default_fin_format  = '%03d_*.dat';
default_fout_format = 'matrix_data_%03d-%03d.mat';

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});

% reset defaults based on config entries
if isfield(config, 'data_directory'); default_data_directory = config.data_directory; end

% parsed arguments override config fields
addParameter(parser, 'cmap', default_cmap);
addParameter(parser, 'marker', default_marker);
addParameter(parser, 'one_panel', default_one_panel);
addParameter(parser, 'scatter_only', default_scatter_only);
addParameter(parser, 'draggable_line', default_draggable_line);
addParameter(parser, 'data_directory', default_data_directory); % can override
addParameter(parser, 'fin_format', default_fin_format); % can override
addParameter(parser, 'fout_format', default_fout_format); % can override
parse(parser, varargin{:});
data_directory = parser.Results.data_directory;
fin_format = parser.Results.fin_format;
fout_format = parser.Results.fout_format;
cmap = parser.Results.cmap;
marker = parser.Results.marker;
one_panel = parser.Results.one_panel;
scatter_only = parser.Results.scatter_only;
draggable_line = parser.Results.draggable_line;

% validate required config fields
required_fields = {'Vbg_col', 'Vtg_col'};
for field = required_fields
    if ~isfield(config, field)
        error('determine_slope requires <%s> in supplied config', char(field));
    end
end

% check for matrix data (loads faster)
fout_name = sprintf(fout_format, fnums(1), fnums(end));
if exist(fout_name, 'file') == 2
    load_raw_data = false;
else
    load_raw_data = true;
end

if load_raw_data
    args = {'data_directory', data_directory, 'fname_format', fin_format};
    DATA = [];
    EE = [];
    NN = [];
    VBG = [];
    VTG = [];
%     xline(ax, 0, 'k--', 'linewidth', 1.5);
    for fnum = fnums
        try
            fprintf('loading %g...\n', fnum);
            d = readcol(fnum, [data_col, config.Vbg_col, config.Vtg_col], args{:});
        catch
            warning('Error reading file %g', fnum);
            continue
        end
        data = d(:,1);
        Vbg = d(:,2);
        Vtg = d(:,3);
        
        % convert to E,n
        f = volts2E_n(Vbg, Vtg, config);
        E = f(:,1); %[V/nm]
        n = f(:,2); %[nm^-2]

        % add data to list
        DATA = [DATA; data];
        EE = [EE; E];
        NN = [NN; n];
        VBG = [VBG; Vbg];
        VTG = [VTG; Vtg];
    end
    
    E_lims = [min(EE), max(EE)];
    n_lims = [min(NN), max(NN)];
    Vbg_lims = [min(VBG), max(VBG)];
    Vtg_lims = [min(VTG), max(VTG)];
    rows = size(d,1);
    cols = length(fnums);

    if ~scatter_only && ~one_panel
        try
            % reshape as matrix
            data = reshape(DATA, rows, [])';
            Vtg = reshape(VTG, rows, [])';
            Vbg = reshape(VBG, rows, [])';
            fprintf('reshaped data\n');
        catch
            % interpolate
            fprintf('interpolating...');
    %         [n,E] = meshgrid(linspace(n_lims(1),n_lims(2),rows), ...
    %                          linspace(E_lims(1),E_lims(2),cols));
    %         data = griddata(EE,NN,DATA,E,n);

            [Vbg,Vtg] = meshgrid(linspace(Vbg_lims(1),Vbg_lims(2),rows), ...
                                 linspace(Vtg_lims(1),Vtg_lims(2),cols));
            data = griddata(VTG,VBG,DATA,Vtg,Vbg);

            fprintf(' done\n');
        end
    end
    
    % save collated data
%     vars = {'DATA', 'EE', 'NN', 'VBG', 'VTG'};
%     save(fout_name, vars{:});
else
    error('not implemented');
    load(fname);
end

% create figure
figure(201);
clf;
if ~one_panel
    ax1 = subplot(1,2,1); hold on;
    ax2 = subplot(1,2,2); hold on;
else
    ax2 = subplot(1,1,1); hold on;
end
fprintf('drawing...\n');

if ~one_panel
    if scatter_only
        scatter(ax1, VBG, VTG, marker, DATA, 's', 'filled');
    else
        im1 = imagesc(ax1, Vbg_lims, Vtg_lims, data);
        set(im1, 'alphadata', ~isnan(data));
    end

    xlabel(ax1, 'V_{bg} (V)');
    ylabel(ax1, 'V_{tg} (V)');
    axis(ax1, 'tight');
    colormap(ax1, cmap);
end

scatter(ax2, NN, EE, marker, DATA, 's', 'filled');
% im2 = imagesc(ax2, minmax(n), minmax(E), data);

xlabel(ax2, 'Density (nm^{-2})');
ylabel(ax2, 'Electric field (V/nm)');
axis(ax2, 'tight');
grid(ax2, 'on');
colormap(ax2, cmap);

if draggable_line
    % add vertical line for slope alignment
    % line = xline(0, 'm--'); % not draggable
    line = plot([0,0], ylim(ax2), 'm--');
    try draggable(line); end
end

% print slope
fprintf('slope = -C_bg / C_tg = -d_tg / d_bg = -%g\n', config.dtg/config.dbg);