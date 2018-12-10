function transistorON()
    smset('Vdd', 0.2, 1);
    smset('Vg', -0.175, 1);
    fprintf('transistor ON\namplifier --> %.4g V\n', cell2mat(smget('Vd')));
end