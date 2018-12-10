function scan_2D(fnum, froot, p0, px, py, Nx, Ny, V1col, V2col, config, varargin)
% perform a 2D parameter scan by looping cosweep(...) with specified start
% and end points along both axes of 2D scan
% this function written by Sergio de la Barrera on Nov 11, 2018
%    p0         <[V1, V2] coordinates of first point in 2D scan>
%    px         <[V1, V2] coordinates of final point along fast axis (plots along x-axis)>
%    py         <[V1, V2] coordinates of final point along slow axis (plots along y-axis)>
%    Nx         <number of points along fast/x-axis>
%    Ny         <number of points along slow/y-axis>
%    V1col		<column used to set V1 values>
%    V2col		<column used to set V2 values>
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and some optionals which can be overridden by varargs:
%                   interval                (see below)
%                   plot_fields             (see below)
%                   limit_condition         (see below)
%
% NO AUTO RAMPING
% PARAMETERS ARE SET INSTANTANEOUSLY ON EACH LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
Vfastrate = 1;
settle_time = 9;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args passed along
validLimitCondition = @(x) validateattributes(x, {'numeric'}, {'numel', 2});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run'})));
addParameter(parser, 'scan_style', 'typewriter', @(x) any(validatestring(x, {'typewriter', 'raster', 'hysteresis'})));
addParameter(parser, 'V1_limits', [-inf,inf], validLimitCondition); % box limits on each separate voltage
addParameter(parser, 'V2_limits', [-inf,inf], validLimitCondition); % (absolute limits beyond which the voltage will not go)

parse(parser, varargin{:});
dry_run     = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run; fprintf('DRY-RUN ONLY\n'); end
raster      = strcmp(parser.Results.scan_style, 'raster');
hysteresis  = strcmp(parser.Results.scan_style, 'hysteresis');
V1_limits   = parser.Results.V1_limits;
V2_limits   = parser.Results.V2_limits;

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
    
    % trim points outside of box limits
    box = V1>=V1_limits(1) & V1<=V1_limits(2) & V2>=V2_limits(1) & V2<=V2_limits(2);
    V1 = V1(box);
    V2 = V2(box);
    
    if isempty(V1) % skip sweeps that fall entirely outside box limits
        fprintf('skipping empty sweep %g/%g\n', ny, Ny);
    else        
        if dry_run
            % record sweep points for plotting
            if ny == 1
                V1_full_scan = V1;
                V2_full_scan = V2;
            else
                V1_full_scan = [V1_full_scan, V1];
                V2_full_scan = [V2_full_scan, V2];
            end
            if hysteresis
                V1_full_scan = [V1_full_scan, flip(V1)];
                V2_full_scan = [V2_full_scan, flip(V2)];
            end
        else
            % fast ramp to start and cosweep along fast axis
            fprintf('ramping...');
            smset({config.channels{V1col}, config.channels{V2col}}, [V1(1), V2(1)], Vfastrate);
            fprintf('settling...\n');
            pause(settle_time);
            cosweep(fnum, froot, V1, V2, V1col, V2col, config, varargin{:});
            if hysteresis
                pause(settle_time);
                fnum = fnum + 1;
                cosweep(fnum, froot, flip(V1), flip(V2), V1col, V2col, config, varargin{:});
            end
        end
        fprintf('completed %g/%g sweeps\n', ny, Ny);
        fnum = fnum + 1;
    end
end

if dry_run
    % simply plot scan itinerary
    figure();
    plot(V1_full_scan, V2_full_scan, '-x');
    xlabel(config.columns{V1col});
    ylabel(config.columns{V2col});
    fprintf('*** DRY-RUN ONLY\n');
else
    % exit
    fprintf('*** %s\t2D scan complete\n', datestr(clock, 'mmm dd HH:MMPM'));
end