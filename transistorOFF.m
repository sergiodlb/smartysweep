function transistorOFF()
    smset({'K2400A.V', 'K2400B.V'}, 0, 1);
    fprintf('transistor OFF\nK2000 --> %.4g V\n', cell2mat(smget('K2000.V')));
end