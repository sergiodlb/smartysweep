function out = balance_capacitance_bridge(config, varargin)
% balances capacitance bridge for subsequent capacitance measurement
% based on capbridgeBen written by Nick Goman and Ben Hunt
% this function written by Sergio de la Barrera on Aug 29, 2017
%    config     structure containing:
%                   channels = {...} (like data_fields)
%                   columns = {...} 
%               and channels specific to capacitance:
%                   Vex_amplitude_channel   <sm channel of excitation>
%                   Vstd_amplitude_channel  <sm channel of standard capacitor>
%                   Vstd_phase_channel      <sm channel of Vstd phase relative to Vex>
%                   time_constant_channel   <sm channel of AC time constant>
%                   Xcol       <column # of off-balance X measurement>
%                   Ycol       <column # of off-balance Y measurement>
%               additional fields required if called by measurement script:
%                   Ccol       <column # to store capacitance measurement>
%                   Lcol       <column # to store capacitance loss>
%               and some optionals which can be overridden by varargs:
%                   Vex                     (see below)
%                   Vstd_range              (see below)
%                   Cstd                    (see below)
% ---- optional parameters (will override duplicate entries in config) ----
%    Vex        <excitation voltage to use; default USE PRESENT VALUE>
%    Vstd_range <range of voltage to use for standard capacitor; default USE PRESENT VALUE>
%    Cstd       <standard capacitance; default = 1>
%    quiet      <BOOL to block text output to stdout; default = false>
%    return_values  <toggles measurement mode (default = false), which will return column numbers and C, Closs values if called by another measurement function>
%                   e.g.
%                     if return_values
%                         out.columns = [config.Ccol, config.Lcol]; % tell parent function where to put data values
%                         out.values  = [Cex, Closs]; % return data values to be stored in the specified columns
%                     else
%                         out = [balance_matrix, Vc0Vex, Vr0Vex, Cex]; % simply return balance_matrix and other measured values
%                     end
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
% 2018-07-20    - major change to output format; if optional parameter
%                 return_values == true, then returned object is a
%                 structure with named fields:
%                    out.columns = [config.Ccol, config.Lcol]
%                    out.values  = [Cex, Closs]
%                 otherwise output structure contains:
%                    out.balance_matrix = balance_matrix; 
%                    out.Vc0Vex         = Vc0Vex;
%                    out.Vr0Vex         = Vr0Vex;
%                    out.Cex            = Cex;
%                    out.Closs          = Closs;
%                    out.error_x        = error_x;
%                    out.error_y        = error_y;
% 2018-07-25    - added text output and associated 'quiet' option
% 2019-04-12    - significant revisions which move the code to handle ZI
%                 behavior with regard to setting output voltages fully
%                 into the ZI driver file
%               - this function (and other related ones) now considers all
%                 AC voltages to be RMS voltages, not peak-to-peak, as
%                 was the case for outputs (but not inputs) before
%               - this also allows interoperability with SR860 lock-ins,
%                 which use RMS voltages for all inputs and outputs
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters that change
% vc = 0.01;	% fraction of output voltage scale
% dvc = 0.49;	% large off-balance variation in fraction of output voltage scale
% vr = 0.01;	% fraction of output voltage scale
% dvr = 0.49;	% large off-balance variation in fraction of output voltage scale
% vc = 0.01;	% fraction of output voltage scale
% dvc = 0.19;	% large off-balance variation in fraction of output voltage scale
% vr = 0.01;	% fraction of output voltage scale
% dvr = 0.19;	% large off-balance variation in fraction of output voltage scale
vc = 0.01;	% fraction of output voltage scale
dvc = 0.94;	% large off-balance variation in fraction of output voltage scale
vr = 0.01;	% fraction of output voltage scale
dvr = 0.94;	% large off-balance variation in fraction of output voltage scale
time_constant_mult  = 10; % how long to wait after changing voltage
default_Vex         = []; % use present value if empty
default_Vstd_range  = []; % use present value if empty
default_Vex_scale   = 1; % external gain (>1) or attenuation (<1)
default_Vstd_scale  = 1; % external gain (>1) or attenuation (<1)
default_Cstd        = 1; % if standard capacitance is unknown
default_offbal_samples = 100; % for averaging (separate from sweep sampling)
default_balance     = true; % wait for true balance at end
default_quiet       = false; % block all text output (other than errors) if true
default_plotting    = false; % plot real-time off-balance voltages while balancing

% validate required config fields
required_fields = {'Vex_amplitude_channel',  'time_constant_channel', ...
                   'Vstd_amplitude_channel', 'Vstd_phase_channel', ...
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
if isfield(config, 'Vex_scale'); default_Vex_scale = config.Vex_scale; end
if isfield(config, 'Vstd_scale'); default_Vstd_scale = config.Vstd_scale; end
if isfield(config, 'Cstd'); default_Cstd = config.Cstd; end

% parsed arguments override config fields
addParameter(parser, 'Vex', default_Vex, validScalarNonNeg); % can override
addParameter(parser, 'Vstd_range', default_Vstd_range, validScalarPos); % can override
addParameter(parser, 'Vex_scale', default_Vex_scale, validScalarNonNeg); % can override
addParameter(parser, 'Vstd_scale', default_Vstd_scale, validScalarNonNeg); % can override
addParameter(parser, 'Cstd', default_Cstd, validScalarPos); % can override
addParameter(parser, 'offbal_samples', default_offbal_samples, validScalarNonNeg);
addParameter(parser, 'balance', default_balance);
addParameter(parser, 'quiet', default_quiet);
addParameter(parser, 'plot', default_plotting);
addParameter(parser, 'return_values', false); % true if called by a measurement script (e.g. to perform rebalanced measurement)

parse(parser, varargin{:});
Vex         = parser.Results.Vex;
Vstd_range  = parser.Results.Vstd_range;
Vex_scale   = parser.Results.Vex_scale;
Vstd_scale  = parser.Results.Vstd_scale;
Cstd        = parser.Results.Cstd;
samples     = parser.Results.offbal_samples;
balance     = parser.Results.balance;
quiet       = parser.Results.quiet;
plotting    = parser.Results.plot;

% check for additional arguments needed if returning data values
return_values = parser.Results.return_values;
if return_values
    required_fields = {'Ccol', 'Lcol'};
    for field = required_fields
        if ~isfield(config, field)
            error('balance_capacitance_bridge_vrms requires <%s> in supplied config in order to return data for measurement', char(field));
        end
    end
end    

if ~isempty(Vex) % user supplied Vex to SET   
    if Vex > 2     % double check with user if input Vex > 2 V
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
    
    % set Vex (rms voltage)
    smset(config.Vex_amplitude_channel, Vex/Vex_scale);%, 0.5);
else % read present values and use those instead
    Vex = cell2mat(smget(config.Vex_amplitude_channel));
end

if ~isempty(Vstd_range) % user supplied Vstd_range to SET
    V_hi = max(sqrt(vr^2 + (vc + dvc)^2), sqrt(vc^2 + (vr + dvr)^2)) * Vstd_range; % highest standard voltage that will be applied
    if V_hi > 2     % double check with user if any V_std > 2 V
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
else
    error('balance_capacitance_bridge_vrms requires <Vstd_range> in supplied config or as optional argument');
end

% measure off-balance voltage components at three points (using fraction of range)
Vcs = [vc, vc, vc + dvc]*Vstd_range; % rms voltages
Vrs = [vr, vr + dvr, vr]*Vstd_range;
L = zeros(2, 3, samples);
if ~quiet; fprintf('balancing'); end
for n = 1:3
    R = sqrt(Vcs(n)^2 + Vrs(n)^2);
    phase = 180 - atan2d(Vrs(n), Vcs(n)); % 4-quadrant tangent function
    smset(config.Vstd_amplitude_channel, R/Vstd_scale);
    smset(config.Vstd_phase_channel, phase);
    pause(time_constant_mult*time_const);
    
    for m = 1:samples
%         while cell2mat(smget('Vd')) > 60e-3 || cell2mat(smget('Vd')) < 50e-3
%             pause(0.1)
%         end
        L(1, n, m) = cell2mat(smget(config.channels{Xcol})); % rms voltages
        L(2, n, m) = cell2mat(smget(config.channels{Ycol}));
        if plotting
            if n == 1 && m == 1
                figure();
                L1 = zeros(1,3*samples);
                L2 = zeros(1,3*samples);
                k = 1;
                L1(k) = L(1,1,1);
                L2(k) = L(2,1,1);
                ax1 = plot(L1(1),'b'); hold all;
                ax2 = plot(L2(1),'r');
            else
                k = k + 1;
                L1(k) = L(1,n,m);
                L2(k) = L(2,n,m);
                set(ax1, 'YData', L1(1:k));
                set(ax2, 'YData', L2(1:k));
                drawnow;
            end
        end
    end
    if ~quiet; fprintf('.'); end
end
L = mean(L, 3); % average over samples
% if ~quiet; disp(L); end;

% convert remaining fractional voltages to real voltage units
Vr = vr*Vstd_range;
Vc = vc*Vstd_range;
dVr = dvr*Vstd_range;
dVc = dvc*Vstd_range;

% the algorithmic part; see Ashoori thesis
Kr1 = (L(1,2) - L(1,1)) / dVr; % real voltage units (but K's are dimensionless)
Kc1 = (L(1,3) - L(1,1)) / dVc;
Kr2 = (L(2,2) - L(2,1)) / dVr;
Kc2 = (L(2,3) - L(2,1)) / dVc;
P = (1 - (Kc1 * Kr2) / (Kr1 * Kc2))^(-1);
Vr0 = Vr + (P / Kr1) * ((Kc1 / Kc2) * L(2,1) - L(1,1)); % all rms voltages
Vc0 = Vc + (P / Kc2) * ((Kr2 / Kr1) * L(1,1) - L(2,1));

% calculate device capacitance (rms voltages)
Cex = Cstd * Vc0 / Vex; % edit on 1/8/2019 for tuning antisymmetric capacitance (can be negative)
Closs = Cstd * Vr0 / Vex;

% % phase-shift corrected version:
% eid = exp(1j * atan2(Vr0, Vc0))
% Cex = Cstd * (Vc0 + Vr0*imag(eid)) / Vex;
% Closs = Cstd * Vr0 * real(eid) / Vex;

% output Vc0/Vex and Vr0/Vex
Vc0Vex = Vc0/Vex;
Vr0Vex = Vr0/Vex;

% set excitation to balance point
R = sqrt(Vc0^2 + Vr0^2);
if R > Vstd_range
    cprintf('red', 'BALANCE POINT Vstd LARGER THAN AVAILABLE RANGE--- INCREASE RANGE OR DROP Vex\n');
    in_range = false;
else
    in_range = true;
    phase = 180 - atan2d(Vr0, Vc0); % 4-quadrant tangent function
    smset(config.Vstd_amplitude_channel, R/Vstd_scale); % go to balance point (if within acceptable Vstd_range)
    smset(config.Vstd_phase_channel, phase);
end

if balance
    % this settle time allows a proper off-balance voltage measurement (ideally zero) following the balancing process
    pause(time_constant_mult*time_const);
    if ~quiet
        if in_range
            fprintf('balanced');
        else
            cprintf('unbalanced!');
        end
    end
    error_x = cell2mat(smget(config.channels{Xcol})); % it also allows determination of the error wrt true balance
    error_y = cell2mat(smget(config.channels{Ycol}));
end
if ~quiet; fprintf('\n'); end

% output result (includes real voltages)
balance_matrix = [Kc1, Kc2, Kr1, Kr2, Vc0, Vr0];

% handle data output for rebalanced measurement versus direct call by user
if return_values
    out.columns = [config.Ccol, config.Lcol]; % tell parent function where to put data values
    out.values  = [Cex, Closs]; % return data values to be stored in the specified columns
else % simply return balance_matrix and other measured values in a structure
    out.balance_matrix = balance_matrix; 
    out.Vc0Vex = Vc0Vex;
    out.Vr0Vex = Vr0Vex;
    out.Cex = Cex;
    out.Closs = Closs;
    out.R = R;
    out.phase = phase;
    if balance % only report error if allowed to balance
        out.error_x = error_x;
        out.error_y = error_y;
    end
end
return