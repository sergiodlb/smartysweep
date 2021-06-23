function out = check_transistor_setpoint(config, varargin)
% attempts to adjust transistor Vg to maintain Vd setpoint value
% works best in linear regime (Vg close to Vg*)
% written by Sergio de la Barrera on 4/10/2019

% defaults
% Vdd_channel = 'Vdd';
Vg_channel  = 'Vg';
Vd_channel  = 'Vd';
tolerance   = 0.02;
delta       = 2; % multiple of Vd-deviation (%) over which to sweep Vg
N           = 21; % number of Vg points to sweep during adjustment
Vg_limit    = -0.5; % lower limit on Vg (check FHX35X datasheet)

% validate required config fields
required_fields = {'Vg', 'Vd_setpoint'};
for field = required_fields
    if ~isfield(config, field)
        error('check_transistor_setpoint requires <%s> in supplied config', char(field));
    end
end

% reset defaults based on config entries
% if isfield(config, 'Vdd_channel'); Vdd_channel = config.Vdd_channel; end
if isfield(config, 'Vg_channel');  Vg_channel  = config.Vg_channel;  end
if isfield(config, 'Vd_channel');  Vd_channel  = config.Vd_channel;  end

% check transistor output
Vd = cell2mat(smget(Vd_channel));
Vd_diff = Vd-config.Vd_setpoint;
if abs(Vd_diff)/config.Vd_setpoint > tolerance
    % sweep a small range of Vg and use the measured Vd values to set a new Vg
    Vd_new = zeros(1, N);
    Vg0    = cell2mat(smget(Vg_channel)); % may have changed from config.Vg due to prior adjustment
    Vgf    = Vg0*(1 - delta*Vd_diff/config.Vd_setpoint); % direction of sweep depends on sign of Vd_diff
    n      = 1;
    Vg_new = linspace(Vg0, Vgf, N);
    for Vg = Vg_new
        smset(Vg_channel, Vg, 1);
        Vd_new(n) = cell2mat(smget(Vd_channel));
        n = n + 1;
    end
    
%     figure();
%     plot(Vg_new, Vd_new);
%     hold on;
    
    % recover set point by linear extrapolation
    p = polyfit(Vg_new, Vd_new, 1);
    Vg_new = (config.Vd_setpoint - p(2))/p(1);
        
%     plot(Vg_new, polyval(p, Vg_new), 's');
    
    % test if within limits
    if Vg_new < 0 && Vg_new > Vg_limit
        smset(Vg_channel, Vg_new, 1);
        
%         pause(0.1);
        fprintf('Vg adjusted to %.4g V,\tamplifier --> %.4g V\n', Vg_new, cell2mat(smget(Vd_channel)));
    else
        warning('Transistor is deviating from setpoint but extrapolated Vg is outside of specified limits --> %.4g V\n', Vg_new);
    end
    out = Vg_new;
    
%     % find optimal point by fitting 3rd-degree polynomial
%     p = polyfit(Vg_new, Vd_new, 3) % f(x) = a*x^3 + b*x^2 + c*x + d --> f'(x) = 3*a*x^2 + 2*b*x + c
%     figure();
%     plot(Vg_new, Vd_new);
%     hold on;
%     plot(Vg_new, polyval(p, Vg_new), 'x');
%     Vg_new = -p(2)/(3*p(1)) % x0 = -2*b/(2*3*a) = -b/(3*a)
%     plot(Vg_new, polyval(p, Vg_new), 's');
else
    out = true;
end
return
