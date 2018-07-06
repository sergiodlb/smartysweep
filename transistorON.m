function transistorON()
    smset('K2400A.V', 0.2, 1);
    smset('K2400B.V', -0.175, 1);
    fprintf('transistor ON\nK2000 --> %.4g V\n', cell2mat(smget('K2000.V')));
end