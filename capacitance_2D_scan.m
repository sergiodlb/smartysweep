function capacitance_2D_scan(fnum, froot, p0, px, py, interval, Nx, Ny, V1col, V2col, Ccol, Lcol, Xcol, Ycol, ex_amplitude_channel, std_amplitude_channel, Vrange_std, std_phase_channel, balance_matrix, data_fields, plot_fields, varargin)
% perform a 2D parameter scan while measuring capacitance
% uses a looping capacitance_cosweep with arbitrary start and end points
% this function written by Sergio de la Barrera on Apr 16, 2018
%    p0         <[V1, V2] coordinates of first point in 2D scan>
%    px         <[V1, V2] coordinates of final point along fast axis (plots along x-axis)>
%    py         <[V1, V2] coordinates of final point along slow axis (plots along y-axis)>
%    Nx         <number of points along fast/x-axis>
%    Ny         <number of points along slow/y-axis>
%    V1col		<column used to set V1 values>
%    V2col		<column used to set V2 values>
%    Ccol		<column that will be used to store capacitance>
%    Lcol		<column that will be used to store capacitance loss>
%    Xcol       <column of device X measurement>
%    Ycol       <column of device Y measurement>
%    ex_amplitude_channel	<sm channel string of AC excitation voltage for device under test>
%    balance_matrix	<matrix of cap-bridge balanced parameters: [Kc1, Kc2, Kr1, Kr2, vc0, vr0]>
%    Cstd       (optional: defaults to 1)<standard capacitance>
%
% NO AUTO RAMPING
% PARAMETERS ARE SET INSTANTANEOUSLY ON EACH LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
Vfastrate = 0.5;
settle_time = 0;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args passed along
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run'})));
addParameter(parser, 'scan_style', 'typewriter', @(x) any(validatestring(x, {'typewriter', 'raster'})));

parse(parser, varargin{:});
dry_run = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run; fprintf('DRY-RUN ONLY\n'); end
raster = strcmp(parser.Results.scan_style, 'raster');

% determine start point of each fast sweep
p0prime = [linspace(p0(1), py(1), Ny); linspace(p0(2), py(2), Ny)];

% initiate 2D scan
for ny = 1:Ny
    % determine end point of fast sweep
    pxprime = px + p0prime(:,ny)' - p0;
    V1 = linspace(p0prime(1,ny), pxprime(1), Nx);
    V2 = linspace(p0prime(2,ny), pxprime(2), Nx);
    if raster && mod(ny-1, 2) % change sweep direction every other scan
        V1 = flip(V1);
        V2 = flip(V2);
    end
    
    if dry_run
        % record sweep points for plotting
        if ny == 1
            V1_full_scan = V1;
            V2_full_scan = V2;
        else
            V1_full_scan = [V1_full_scan, V1];
            V2_full_scan = [V2_full_scan, V2];
        end
    else
        % fast ramp to start and cosweep along fast axis
        smset({data_fields.channels{V1col}, data_fields.channels{V2col}}, [V1(1), V2(1)], Vfastrate);
        pause(settle_time);
        capacitance_cosweep(fnum, froot, V1, V2, interval, V1col, V2col, Ccol, Lcol, Xcol, Ycol, ex_amplitude_channel, std_amplitude_channel, Vrange_std, std_phase_channel, balance_matrix, data_fields, plot_fields, varargin{:});
    end
    fprintf('completed %g/%g sweeps\n', ny, Ny);
    fnum = fnum + 1;
end

if dry_run
    % simply plot scan itinerary
    figure();
    plot(V1_full_scan, V2_full_scan, '-x');
    xlabel(data_fields.columns{V1col});
    ylabel(data_fields.columns{V2col});
else
    % exit
    fprintf('*** %s\tcapacitance 2D scan complete\n', datestr(clock, 'mmm dd HH:MMPM'));
end