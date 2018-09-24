%% SM instrument initialization
clear all;
close all;
instrreset;
global smdata;

% GPIB addresses
GPIB_BOARD = 'ni';
BOARD_NUM = 0;

SR860_GPIB = [4,6];
K2400_GPIB = [24];
K2000_GPIB = 16;
K6500_GPIB = 17;
AMI430_COM = 'COM4';
AMI430_BAUD = 115200;
SeekatUNO_COM = 'COM5';
SeekatUNO_BAUD= 9600;
%LS372_GPIB = 12;

labels = 'ABCD';

%% add AMI430 magnet controller
try
    ind = smloadinst('AMI430', [], AMI430_COM, 'BaudRate', AMI430_BAUD);
    
    % open serial communication
    smopen(ind);
    smdata.inst(ind).name = 'AMI430'; 
    
    % add channels
    smaddchannel('AMI430', 'FIELD', 'AMI430.B');
    smaddchannel('AMI430', 'RATE1', 'AMI430.rate1');
    smaddchannel('AMI430', 'RATE2', 'AMI430.rate2');
    smaddchannel('AMI430', 'RAMP', 'AMI430.ramp');
    smaddchannel('AMI430', 'PAUSE', 'AMI430.pause');
    smaddchannel('AMI430', 'PSWCH', 'AMI430.pswitch');
    smaddchannel('AMI430', 'STATE', 'AMI430.state');
    
catch err
    fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n']);
end

%% Add QD 6000 to rack
try
    ind = smloadinst('QD6000',[], GPIB_BOARD, BOARD_NUM, QD6000_GPIB);
    
    smopen(ind);
    smdata.inst(ind).name = 'PPMS'; 
    
    smaddchannel('PPMS', 'TEMP', 'PPMS.T', [0, Inf, Inf, 1]); % [min, max, rate_max, factor]
    smaddchannel('PPMS', 'TRATE', 'PPMS.Trate', [0.01, 20, Inf, 1]);
    smaddchannel('PPMS', 'DFIELD1SHOT', 'PPMS.B');
    smaddchannel('PPMS', 'BRATE', 'PPMS.Brate', [10.3, 190.5, Inf, 1]);
%     smaddchannel('PPMS', 'PFIELD1SHOT', 'PPMS.Pfield'); % persistent field; we don't typically use this
    smaddchannel('PPMS', 'HELIUMLEVEL', 'PPMS.He');
    smaddchannel('PPMS', 'STANDBY', 'PPMS.standby'); % puts PPMS in standby mode (to reduce He boil-off); set-only
    
catch err
   fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n'])
end

%% add SR860 lock-ins
n = 0;
for GPIB = SR860_GPIB
    try
        n = n + 1;
        ind = smloadinst('SR860', [], GPIB_BOARD, BOARD_NUM, GPIB);

        % open GPIB communication
        smopen(ind);
        if length(SR860_GPIB) > 1
            name = ['SR860', labels(n)];
        else
            name = 'SR860';
        end
        smdata.inst(ind).name = name;

        % add channels
        smaddchannel(name, 'X', [name, '.X']); % fixed spaces (unicode 32) issue by replacement (unicode 0) --> inst.channels(uint8(inst.channels) == 32) = char(0)
        smaddchannel(name, 'Y', [name, '.Y']);
        smaddchannel(name, 'R', [name, '.R']);
        smaddchannel(name, 'THETA', [name, '.phase']);
        smaddchannel(name, 'FREQ', [name, '.freq']);
        smaddchannel(name, 'VREF', [name, '.V']);
        smaddchannel(name, 'SOFF', [name, '.Vdc']);
        smaddchannel(name, 'SCAL', [name, '.sensitivity']);
        smaddchannel(name, 'TAU', [name, '.timeconstant']);

    catch err
       fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n']);
    end
end

%% Add ZI HF2LI lock-in to rack
try
    ind = smloadinst('ZIHF2LI');
    ziDAQ('connect');
    smdata.inst(ind).name = 'ZI'; 
    
    % add channels
    smaddchannel('ZI', 'DEMOD1SAMPX', 'ZI.1X');
    smaddchannel('ZI', 'DEMOD1SAMPY', 'ZI.1Y');
    smaddchannel('ZI', 'DEMOD1PHASE', 'ZI.1phase');
    smaddchannel('ZI', 'DEMOD1TC', 'ZI.1timeConstant');
%     smaddchannel('ZI', 'IN1RANGE', 'Input 1 Range');
    
    smaddchannel('ZI', 'DEMOD2SAMPX', 'ZI.2X');
    smaddchannel('ZI', 'DEMOD2SAMPY', 'ZI.2Y');
    smaddchannel('ZI', 'DEMOD2PHASE', 'ZI.2phase');
    smaddchannel('ZI', 'DEMOD2TC', 'ZI.2timeConstant');

    smaddchannel('ZI', 'OSC1F', 'ZI.1freq');
    smaddchannel('ZI', 'OUT1AMP1', 'ZI.OUT1AMP1');
%     smaddchannel('ZI', 'OUT1AMP1ON', 'Output 1 amp 1 on', [0, 1, Inf, 1]);    
%     smaddchannel('ZI', 'OUT1ON', 'Output 1 on', [0, 1, Inf, 1]);
    smaddchannel('ZI', 'OUT1RANGE', 'ZI.OUT1RANGE');
%     
     smaddchannel('ZI', 'OSC2F', 'ZI.2freq');
     smaddchannel('ZI', 'OUT2AMP2', 'ZI.OUT2AMP2');
%     smaddchannel('ZI', 'OUT2AMP2ON', 'Output 2 amp 2 on', [0, 1, Inf, 1]);
%     smaddchannel('ZI', 'OUT2ON', 'Output 2 on', [0, 1, Inf, 1]);
    smaddchannel('ZI', 'OUT2RANGE', 'ZI.OUT2RANGE');
    
%     smaddchannel('ZI', 'BALBRIDGE', 'Balance bridge', [0, 1, Inf, 1]);
%     smaddchannel('ZI', 'DELTAC', 'Capacitance');
    
catch err
    fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n']);
end

%% add SeekatUNO 8-channel DC voltage source
try
    ind = smloadinst('SeekatUNO', [], SeekatUNO_COM, 'BaudRate', SeekatUNO_BAUD);

    % open communication
    smopen(ind) 
    smdata.inst(ind).name = 'SeekatUNO'; 
    
    % add channels
    smaddchannel('SeekatUNO', 'V1', 'Vdd');
    smaddchannel('SeekatUNO', 'V2', 'Vg');
    smaddchannel('SeekatUNO', 'V3', 'SeekatUNO.V3');
    smaddchannel('SeekatUNO', 'V4', 'SeekatUNO.V4');
    smaddchannel('SeekatUNO', 'V5', 'SeekatUNO.V5');
    smaddchannel('SeekatUNO', 'V6', 'SeekatUNO.V6');
    smaddchannel('SeekatUNO', 'V7', 'SeekatUNO.V7');
    smaddchannel('SeekatUNO', 'V8', 'SeekatUNO.V8'); % limits appear to be +/- 10 V
    
    smset('Vdd', 0);
    smset('Vg', 0);
    smset('SeekatUNO.V3', 0);
    smset('SeekatUNO.V4', 0);
    smset('SeekatUNO.V5', 0);
    smset('SeekatUNO.V6', 0);
    smset('SeekatUNO.V7', 0);
    smset('SeekatUNO.V8', 0);
catch err 
   fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n'])
end

%% loading K2400s
n = 0;
for GPIB = K2400_GPIB
    try
        n = n + 1;
        ind = smloadinst('Keithley2400', [], GPIB_BOARD, BOARD_NUM, GPIB);
        
        % open GPIB communicationZI
        smopen(ind);
        if length(K2400_GPIB) > 1
            name = ['K2400', labels(n)];
        else
            name = 'K2400';
        end
        smdata.inst(ind).name = name;

        % add channels
        smaddchannel(name, 'VOLT', [name, '.V']);
        smaddchannel(name, 'CURRENT', [name, '.I']);
        smaddchannel(name, 'COMPLIANCE', [name, '.compliance']);
        smaddchannel(name, 'ISSOURCEVOLT', [name, '.sourceVolt']);
        smaddchannel(name, 'OUTPUTON', [name, '.outputOn']);
        smaddchannel(name, 'RESISTANCE', [name, '.R']);
        smaddchannel(name, 'VRANGE', [name, '.range']);
        
        % ensure voltage sourcing mode, turn output ON
        if ~cell2mat(smget([name, '.sourceVolt']))
            smset([name, '.sourceVolt'], 1);
        end
        if ~cell2mat(smget([name, '.outputOn']))
            smset([name, '.outputOn'], 1);
        end
    
    catch err
       fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n']);
    end
end

%% add Keithley 2000 digital multimeter (using Keithley 2002 driver)
try
    ind = smloadinst('Keithley2002', [], GPIB_BOARD, BOARD_NUM, K2000_GPIB);
     
    % open GPIB communication
    smopen(ind) 
    smdata.inst(ind).name = 'K2000'; 
    
    % add channels
    smaddchannel('K2000', 'VOLT', 'K2000.V', [-Inf, Inf, Inf, 1]);
    
catch err
   fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n'])
end

%% add Keithley DMM6500 digital multimeter
try
    ind = smloadinst('Keithley6500', [], GPIB_BOARD, BOARD_NUM, K6500_GPIB);
     
    % open GPIB communication
    smopen(ind) 
    smdata.inst(ind).name = 'K6500'; 
    
    % add channels
    smaddchannel('K6500', 'VOLT', 'K6500.V', [-Inf, Inf, Inf, 1]);
    smaddchannel('K6500', 'CURRENT', 'K6500.I', [-Inf, Inf, Inf, 1]);
    
catch err
   fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n'])
end

%% add LakeShore 372
try
    ind = smloadinst('LS372', [], GPIB_BOARD, BOARD_NUM, LS372_GPIB);
     
    % open GPIB communication
    smopen(ind) 
    smdata.inst(ind).name = 'LS372'; 
    
    % add channels
    smaddchannel('LS372', 'T1', 'LS372.T1');
    smaddchannel('LS372', 'T2', 'LS372.T2');
    smaddchannel('LS372', 'T3', 'LS372.T3');
    smaddchannel('LS372', 'T5', 'LS372.T5');
    smaddchannel('LS372', 'T6', 'LS372.T6');
    smaddchannel('LS372', 'T9', 'LS372.T'); % let's just call the probe channel 'T' for simplicity
    smaddchannel('LS372', 'OUTMODE', 'LS372.outmode'); % set 1 for closed loop PID and 0 for OFF
    smaddchannel('LS372', 'RANGE', 'LS372.range'); % heater range: 0 = off, 1 = 31.6 µA, 2 = 100 µA, 3 = 316 µA, 4 = 1.00 mA, 5 = 3.16 mA, 6 = 10.0 mA, 7 = 31.6 mA, 8 = 100 mA
    smaddchannel('LS372', 'SETP', 'LS372.setpoint');
    smaddchannel('LS372', 'RAMP', 'LS372.ramp'); % turn setpoint ramping ON/OFF: 1=ON, 0=OFF
    smaddchannel('LS372', 'RAMPRATE', 'LS372.rate'); % set/get setpoint ramp rate
    smaddchannel('LS372', 'AUTOSCAN', 'LS372.autoscan'); % turn autoscanning ON/OFF: 1=ON, 0=OFF
    
catch err
   fprintf(['*ERROR* problem with connecting to device\n' err.identifier ': ' err.message '\n'])
end