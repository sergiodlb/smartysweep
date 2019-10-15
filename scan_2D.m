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
Vfastrate = 0.2;
settle_time = 3;

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args passed along
validLimitCondition = @(x) validateattributes(x, {'numeric'}, {'numel', 2});
validFunction = @(x) validateattributes(x, {'function_handle'}, {});
addOptional(parser, 'dry_run', false, @(x) any(validatestring(x, {'dry_run'})));
addParameter(parser, 'scan_style', 'typewriter', @(x) any(validatestring(x, {'typewriter', 'raster', 'hysteresis'})));
addParameter(parser, 'V1_limits', [-inf,inf], validLimitCondition); % box limits on each separate voltage
addParameter(parser, 'V2_limits', [-inf,inf], validLimitCondition); % (absolute limits beyond which the voltage will not go)
addParameter(parser, 'mask_function', [], validFunction); % function which takes f(V1, V2, config) and returns BOOLEAN ARRAY indicating valid voltage values for sweeping

parse(parser, varargin{:});
dry_run     = parser.Results.dry_run; % will behave as true if dry_run is string
if dry_run; fprintf('DRY-RUN ONLY\n'); end
raster      = strcmp(parser.Results.scan_style, 'raster');
hysteresis  = strcmp(parser.Results.scan_style, 'hysteresis');
V1_limits   = parser.Results.V1_limits;
V2_limits   = parser.Results.V2_limits;
mask_function = parser.Results.mask_function;

% determine start point of each fast sweep
p0prime = [linspace(p0(1), py(1), Ny); linspace(p0(2), py(2), Ny)];

% set up unequal hysteresis sweeps
if hysteresis && length(Nx) > 1
    Nx_back = Nx(2);
    Nx      = Nx(1);
    unequal_hysteresis = true;
else
    unequal_hysteresis = false;
end

% initiate 2D scan
for ny = 1:Ny
    % determine end point of fast sweep
    pxprime = px + p0prime(:,ny)' - p0;
    V1 = linspace(p0prime(1,ny), pxprime(1), Nx);
    V2 = linspace(p0prime(2,ny), pxprime(2), Nx);
    if raster && mod(ny-1, 2) % change sweep direction every other scan
        V1 = flip(V1);
        V2 = flip(V2);
    elseif hysteresis
        if unequal_hysteresis
            V1_back = linspace(pxprime(1), p0prime(1,ny), Nx_back);
            V2_back = linspace(pxprime(2), p0prime(2,ny), Nx_back);
        else
            V1_back = flip(V1);
            V2_back = flip(V2);
        end
    end
    
    % trim points outside of box limits
    box = V1>=V1_limits(1) & V1<=V1_limits(2) & V2>=V2_limits(1) & V2<=V2_limits(2);
    if ~isempty(mask_function)
        mask = mask_function(V1, V2, config);
        box = box & mask;
    end
    V1 = V1(box);
    V2 = V2(box);
%     if unequal_hysteresis
    if hysteresis
        box = V1_back>=V1_limits(1) & V1_back<=V1_limits(2) & V2_back>=V2_limits(1) & V2_back<=V2_limits(2);
        if ~isempty(mask_function)
            mask = mask_function(V1_back, V2_back, config);
            box = box & mask;
        end
        V1_back = V1_back(box);
        V2_back = V2_back(box);
    end       
    
    if isempty(V1) % skip sweeps that fall entirely outside box limits
        fprintf('skipping empty sweep %g/%g\n', ny, Ny);
    else        
        if dry_run
            % record sweep points for plotting
            if exist('V1_full_scan')
                V1_full_scan = [V1_full_scan, V1];
                V2_full_scan = [V2_full_scan, V2];
            else % first time through, initialize
                V1_full_scan = V1;
                V2_full_scan = V2;
            end
            if hysteresis
                V1_full_scan = [V1_full_scan, V1_back];
                V2_full_scan = [V2_full_scan, V2_back];
            end
        else
            % fast ramp to start and cosweep along fast axis
            fprintf('ramping...');
            smset({config.channels{V1col}, config.channels{V2col}}, [V1(1), V2(1)], Vfastrate);
            fprintf('settling...\n');
            pause(settle_time);
            cosweep_successful = cosweep(fnum, froot, V1, V2, V1col, V2col, config, varargin{:});
            if hysteresis
                if cosweep_successful
                    pause(settle_time);
                else % if, for example, a limit condition was reached in last cosweep, start new cosweep at present values
                    [zero, limit_index] = min(abs(cell2mat(smget(config.channels{V1col}))-V1_back)); % finds closest element to present value
                    V1_back = V1_back(limit_index:end); % truncate both vectors starting at present value
                    V2_back = V2_back(limit_index:end);
                end
                fnum = fnum + 1;
                cosweep(fnum, froot, V1_back, V2_back, V1col, V2col, config, varargin{:});
            end
        end
        fprintf('%s\tcompleted %g/%g sweeps\n', datestr(clock, 'mmm dd HH:MMPM'), ny, Ny);
        fnum = fnum + 1;
    end
end

if dry_run
    % simply plot scan itinerary
    figure();
    l = plot(V1_full_scan, V2_full_scan, '-x');
    hold on;
    plot(V1_full_scan(1), V2_full_scan(1), 'o', 'markerfacecolor', get(l,'color'), 'markeredgecolor', 'w');
    plot(V1_full_scan(end), V2_full_scan(end), 'o', 'markerfacecolor', 'w', 'markeredgecolor', get(l,'color'));
    hold off;
    xlabel(config.columns{V1col});
    ylabel(config.columns{V2col});
    fprintf('*** DRY-RUN ONLY\n');
else
    % exit
    fprintf('*** %s\t2D scan complete\n', datestr(clock, 'mmm dd HH:MMPM'));
end