function transistorOFF()
    smset({'Vdd', 'Vg'}, 0, 1);
    fprintf('transistor OFF\nK2000 --> %.4g V\n', cell2mat(smget('K2000.V')));
end