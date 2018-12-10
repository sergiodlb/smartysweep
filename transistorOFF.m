function transistorOFF()
    smset({'Vdd', 'Vg'}, 0, 1);
    fprintf('transistor OFF\namplifier --> %.4g V\n', cell2mat(smget('Vd')));
end