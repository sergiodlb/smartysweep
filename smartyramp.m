function smartyramp(channels, values, ramprate, varargin)
% ramprate  = 0.2; %[Hz]
ramppause = 100e-3; %[s]
default_quiet = false; % block all text output (other than errors) if true

% deal with optional arguments
parser = inputParser;
parser.KeepUnmatched = true; % other args ignored
addParameter(parser, 'quiet', default_quiet);
parse(parser, varargin{:});
quiet = parser.Results.quiet;

% get present values
values0 = cell2mat(smget(channels));
if length(values0) > length(values)
    if length(values) > 1
        error('Must supply either one set point or an equal number of set points as channels!');
    else
        values = values*ones(size(values0));
    end
end

% determine value differences
dvalues = values - values0;

% determine number of steps
nsteps = round(abs(dvalues)/(ramprate*ramppause));
nsteps = max(nsteps);

setpoints = [];
for n = 1:length(values0)
    setpoints = [setpoints; linspace(values0(n), values(n), nsteps)];
end

if ~isempty(setpoints) && ~quiet
    fprintf('Chan:\t%s\n', sprintf('%-12s', channels{:}));
    fprintf('Init:\t%s\n', sprintf('%-12g', values0'));
end

for setpoint = setpoints
    smset(channels, setpoint);
    if ~quiet
        if exist('line_bytes', 'var')
            fprintf(repmat('\b', 1, line_bytes));
        end
        line_bytes = fprintf('Now:\t%s\n', sprintf('%-12g', setpoint'));
    end
    pause(ramppause);
end