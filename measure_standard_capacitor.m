function [Cstd, freq, Vc0Vex, Vr0Vex] = measure_standard_capacitor(fstart, fend, Npoints, config, varargin)
% measure standard capacitor by sweeping frequency of excitation and
% rebalancing; the standard capacitance is deduced from the slope of the
% resistive component of the output voltage used to balance (linear regime)
%    fstart     <starting frequency for sweep>
%    fend       <ending frequency for sweep>
%    Npoints    <number of frequency points>
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and channels specific to capacitance:
%                   Vex_amplitude_channel   <sm channel of excitation>
%                   Vex_range_channel       <sm channel of Vex range; will use this value to set voltage and range>
%                   Vstd_amplitude_channel  <sm channel of standard capacitor>
%                   Vstd_range_channel      <sm channel of Vstd range; will use to set range and choose off-balance values within this range>
%                   Vstd_phase_channel      <sm channel of Vstd phase relative to Vex>
%                   time_constant_channel   <sm channel of AC time constant>
%                   frequency_channel       <sm channel of AC excitation frequency>
%                   Xcol       <column # of off-balance X measurement>
%                   Ycol       <column # of off-balance Y measurement>
%                   Rchip      <resistance of on-chip high impendance resistor; required argument for measuring standard capacitor>
%               and some optionals which can be overridden by varargs:
%                   Vex                     (see below)
%                   Vstd_range              (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    Vex        <excitation voltage to use; default USE PRESENT VALUE>
%    Vstd_range <range of voltage to use for standard capacitor; default USE PRESENT VALUE>
%    log_scale  <BOOL toggle log scale for frequency sweep>
%
% 2018-04-24    - vast reduction in number of positional arguments; now
%                 using fields in "config" structure to convey capacitance
%                 settings
%               - added optional arguments for basic execution options with
%                 ability to override options in config structure by
%                 choosing name-value pair as optional vararg to function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
default_log_scale   = false; % toggle log x-scale and log spacing of freq points
% time_constant_mult = 10; % how long to wait after setting frequency; EDIT: wait time is built-in to balance_capacitance_bridge()
default_Rchip       = []; % throw error if Rchip is not specified either in config or as optional argument

% validate required config fields (some required by balance_capacitance_bridge)
required_fields = {'time_constant_channel', 'frequency_channel'}; % Rchip handled separately
for field = required_fields
    if ~isfield(config, field)
        error('measure_standard_capacitor requires <%s> in supplied config', char(field));
    end
end

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});
addParameter(parser, 'log_scale', default_log_scale);

if isfield(config, 'Rchip'); default_Rchip = config.Rchip; end
addParameter(parser, 'Rchip', default_Rchip, validScalarPos); % override config.Rchip if given
parse(parser, varargin{:});

% validate on-chip resistance value provided
Rchip = parser.Results.Rchip;
if isempty(Rchip)
    error('measure_standard_capacitor requires Rchip in supplied config or as optional argument');
end

% toggle frequency log scale
log_scale = parser.Results.log_scale;
if log_scale
    freq = logspace(log10(fstart), log10(fend), Npoints);
else
    freq = linspace(fstart, fend, Npoints);
end

% pre-allocate measurement arrays
Vc0Vex = zeros(1, Npoints);
Vr0Vex = zeros(1, Npoints);

% loop over frequencies
for n = 1:Npoints
    smset(config.frequency_channel, freq(n));
    time_constant = 3/freq(n);
    smset(config.time_constant_channel, time_constant); % SET time_constant based on freq
%     pause(time_constant_mult*time_constant); % wait time is built-in to balance_capacitance_bridge()

    [balance_matrix, Vc0Vex(n), Vr0Vex(n), Cex] = balance_capacitance_bridge(config, varargin{:});
    if n == 1
        % create figure first time
        figure();
        ax(1) = plot(freq(1:n), Vc0Vex(1:n), '-x');
        hold all;
        ax(2) = plot(freq(1:n), Vr0Vex(1:n), '-o');
        xlabel('frequency (Hz)');
        ylabel('V_{std}/V_{ex}');
        hold off;
        grid on;
        legend({'V_{c0}', 'V_{r0}'});
        legend('boxoff');
        if log_scale
            set(gca, 'XScale', 'log');
    %         set(gca, 'XScale', 'log', 'YScale', 'log');
        end
    else
        % update existing plots with new data
        try
            set(ax(1), 'XData', freq(1:n), 'YData', Vc0Vex(1:n));
            set(ax(2), 'XData', freq(1:n), 'YData', Vr0Vex(1:n));
            drawnow;
        catch ploterr
            % exit gracefully if plot closed
            if strcmp(ploterr.identifier, 'MATLAB:class:InvalidHandle')
                verb = 'canceled';
                break;
            else
                rethrow(ploterr);
            end
        end
    end
end

% calculate standard capacitance and output result
slope = mean(diff(Vr0Vex)./diff(freq));
fprintf('average slope --> %.4g us\n', slope*1e6);
Cstd = -slope/(2*pi*Rchip);
fprintf('Cstd = slope/(2*pi*Rchip) --> %.4g pF\n', Cstd*1e12);
end