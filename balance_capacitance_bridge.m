function [balance_matrix, Vc0Vex, Vr0Vex, Cex] = balance_capacitance_bridge(config, varargin)
% balances capacitance bridge for subsequent capacitance measurement
% based on capbridgeBen written by Nick Goman and Ben Hunt
% this function written by Sergio de la Barrera on Aug 29, 2017
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
%                   Xcol       <column # of off-balance X measurement>
%                   Ycol       <column # of off-balance Y measurement>
%               and some optionals which can be overridden by varargs:
%                   Vex                     (see below)
%                   Vstd_range              (see below)
%                   Cstd                    (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    Vex        <excitation voltage to use; default USE PRESENT VALUE>
%    Vstd_range <range of voltage to use for standard capacitor; default USE PRESENT VALUE>
%    Cstd       <standard capacitance; default = 1>
%
% 2017-09-08    - heavily revised to use actual voltage units for comparisons
%                 between Vex and std voltages; user provides excitation
%                 voltage (range set based on this value) and RANGE of std
%                 voltages; changed atan --> atan2 (4-quadrant fn) since
%                 former did not appear to balance correctly
% 2018-03-27    - many changes to ensure proper units of each quantity;
%                 several remaining calculations appeared to use fractional 
%                 voltages instead of proper voltage units
% 2018-04-24    - vast reduction in number of positional arguments; now
%                 using fields in "config" structure to convey capacitance
%                 settings
%               - added optional arguments for basic execution options with
%                 ability to override options in config structure by
%                 choosing name-value pair as optional vararg to function
% 2018-07-05    - added Vstd_lockin_range which is separate from Vstd_range
%                 allowing the latter to be less than the available lockin
%                 output range (used during parameter search)
%                 *final balance voltage may be <= Vstd_lockin_range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
vc = 0.01;	% fraction of output voltage scale
dvc = 0.49;	% large off-balance variation in fraction of output voltage scale
vr = 0.01;	% fraction of output voltage scale
dvr = 0.49;	% large off-balance variation in fraction of output voltage scale
samples = 50; % for averaging
time_constant_mult = 10; % how long to wait after changing voltage
default_Vex         = []; % use present value if empty
default_Vstd_range  = []; % use present value if empty
default_Cstd        = 1; % if standard capacitance is unknown

% validate required config fields
required_fields = {'Vex_amplitude_channel',  'Vex_range_channel',  'time_constant_channel', ...
                   'Vstd_amplitude_channel', 'Vstd_range_channel', 'Vstd_phase_channel', ...
                   'Xcol', 'Ycol'};
for field = required_fields
    if ~isfield(config, field)
        error('balance_capacitance_bridge requires <%s> in supplied config', char(field));
    end
end
Xcol = config.Xcol;
Ycol = config.Ycol;
time_const = cell2mat(smget(config.time_constant_channel));

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
validScalarNonNeg = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
validScalarPos = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'});

% reset defaults based on config entries
if isfield(config, 'Vex'); default_Vex = config.Vex; end
if isfield(config, 'Vstd_range'); default_Vstd_range = config.Vstd_range; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end

% parsed arguments override config fields
addParameter(parser, 'Vex', default_Vex, validScalarNonNeg); % can override
addParameter(parser, 'Vstd_range', default_Vstd_range, validScalarPos); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override

parse(parser, varargin{:});
Vex = parser.Results.Vex;
Vstd_range = parser.Results.Vstd_range;
Cstd = parser.Results.Cstd;

if ~isempty(Vex) % user supplied Vex to SET
    % set excitation voltage and range based on Vex input
    if Vex < 10e-3       % use 10 mV range
        Vex_range = 10e-3;
    elseif Vex < 100e-3  % use 100 mV range
        Vex_range = 100e-3;
    elseif Vex < 1 % use 1 V range
        Vex_range = 1;
    else                % use 10 V range
        Vex_range = 10;
        if Vex > 2;     % double check with user if input Vex > 2 V
            msg = sprintf('Selected excitation voltage is %.3g V --> do you want to continue?', Vex);
            disp(msg);
            proceed = menu(msg, 'yes', 'no');
            if proceed == 0 || proceed == 2
                balance_matrix = [];
                Vc0Vex = nan;
                Vr0Vex = nan;
                Cex = nan;
                return
            end
        end
    end

    % set voltage based on computed fraction of range and THEN change range
    vd = Vex/Vex_range;
    if Vex_range > cell2mat(smget(config.Vex_range_channel))
        smset(config.Vex_amplitude_channel, vd);
        smset(config.Vex_range_channel, Vex_range); % avoids jumping to 10 V, for example
    else
        smset(config.Vex_range_channel, Vex_range); % can safely adjust range first
        smset(config.Vex_amplitude_channel, vd); % in case Vd_new > Vd_old
    end
else % read present values and use those instead
    Vex = cell2mat(smget(config.Vex_amplitude_channel));
    Vex_range = cell2mat(smget(config.Vex_range_channel));
    vd = Vex/Vex_range;
end

if ~isempty(Vstd_range) % user supplied Vstd_range to SET
    % set standard voltage and range based on Vstd_range input
%     if ~(Vstd_range == 10^ceil(log10(Vstd_range)))
%         disp('Standard voltage range does not correspond to 10 mV, 100 mV, 1 V, or 10 V --> increasing to next highest range');
%         Vstd_range = 10^ceil(log10(Vstd_range));
%     end
    Vstd_lockin_range = 10^ceil(log10(Vstd_range));
    V_hi = max(sqrt(vr^2 + (vc + dvc)^2), sqrt(vc^2 + (vr + dvr)^2)) * Vstd_range; % highest standard voltage that will be applied
    if V_hi > 2;     % double check with user if any V_std > 2 V
        msg = sprintf('At least one standard voltage will be up to %.3g V --> do you want to continue?', V_hi);
        disp(msg);
        proceed = menu(msg, 'yes', 'no');
        if proceed == 0 || proceed == 2
            balance_matrix = [];
            Vc0Vex = nan;
            Vr0Vex = nan;
            Cex = nan;
            return
        end
    end
    
    if Vstd_lockin_range > cell2mat(smget(config.Vstd_range_channel))
        smset(config.Vstd_amplitude_channel, vc);
        smset(config.Vstd_range_channel, Vstd_lockin_range); % avoids jumping to 10 V, for example
    else
        smset(config.Vstd_range_channel, Vstd_lockin_range); % can safely adjust range first
        smset(config.Vstd_amplitude_channel, vc); % in case Vc_new > Vc_old
    end
else % read present value and use that instead
    Vstd_range = cell2mat(smget(config.Vstd_range_channel));
end

% measure off-balance voltage components at three points (using fraction of range)
vcs = [vc, vc, vc + dvc]*Vstd_range/Vstd_lockin_range;
vrs = [vr, vr + dvr, vr]*Vstd_range/Vstd_lockin_range;
L = zeros(2, 3, samples);
for n = 1:3
    r = sqrt(vcs(n)^2 + vrs(n)^2);
%     phase = 180 - atand(vrs(n) / vcs(n));
    phase = 180 - atan2d(vrs(n), vcs(n)); % 4-quadrant tangent function
    smset(config.Vstd_amplitude_channel, r);
    smset(config.Vstd_phase_channel, phase);
    pause(time_constant_mult*time_const);
    
    for m = 1:samples
        % The measured voltage is 1/1.005 times the output Vrms
        % That is, it's 1 / (1.005 * sqrt(2)) times the output Vpp
        % So the output is really 1.005 * sqrt(2) * the measured input
        L(1, n, m) = 1.005 * sqrt(2) * cell2mat(smget(config.channels{Xcol})); % real voltage units
        L(2, n, m) = 1.005 * sqrt(2) * cell2mat(smget(config.channels{Ycol}));
%         L(1, n, m) = cell2mat(smget(config.channels{Xcol}));
%         L(2, n, m) = cell2mat(smget(config.channels{Ycol}));
    end
end
L = mean(L, 3); % average over samples

% convert remaining fractional voltages to real voltage units
Vr = vr*Vstd_range;
Vc = vc*Vstd_range;
dVr = dvr*Vstd_range;
dVc = dvc*Vstd_range;

% the algorithmic part; see Ashoori thesis
Kr1 = (L(1,2) - L(1,1)) / dVr; % real voltage units
Kc1 = (L(1,3) - L(1,1)) / dVc;
Kr2 = (L(2,2) - L(2,1)) / dVr;
Kc2 = (L(2,3) - L(2,1)) / dVc;
P = (1 - (Kc1 * Kr2) / (Kr1 * Kc2))^(-1);
% vr0 = vr + (P / Kr1) * ((Kc1 / Kc2) * L(2,1) - L(1,1)); % not sure about these units
% vc0 = vc + (P / Kc2) * ((Kr2 / Kr1) * L(1,1) - L(2,1));
Vr0 = Vr + (P / Kr1) * ((Kc1 / Kc2) * L(2,1) - L(1,1)); % should be in proper voltage units
Vc0 = Vc + (P / Kc2) * ((Kr2 / Kr1) * L(1,1) - L(2,1));

% calculate device capacitance (requires real voltages)
% Vc0 = vc0 * Vstd_range;
% Vr0 = vr0 * Vstd_range;
Cex = Cstd * abs(Vc0 / Vex);

% output Vc0/Vex and Vr0/Vex (requires real voltages)
Vc0Vex = Vc0/Vex;
Vr0Vex = Vr0/Vex;

% set excitation to balance point (using fraction of lock-in output range)
vr0 = Vr0/Vstd_lockin_range;
vc0 = Vc0/Vstd_lockin_range;
r = sqrt(vc0^2 + vr0^2);
if r > 1
    cprintf('red', 'BALANCE POINT V_std LARGER THAN AVAILABLE RANGE--- INCREASE RANGE OR DROP Vex\n');
end
% phase = 180 - atand(vr0 / vc0);
phase = 180 - atan2d(vr0, vc0); % 4-quadrant tangent function
smset(config.Vstd_amplitude_channel, r);
smset(config.Vstd_phase_channel, phase); 

% output result (includes real voltages)
balance_matrix = [Kc1, Kc2, Kr1, Kr2, Vc0, Vr0];
end

