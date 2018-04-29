function transistorON()
    smset('Vdd', 0.2, 1);
    smset('Vg', -0.177, 1); % with BG floating
    fprintf('transistor ON\nK2000 --> %.4g V\n', cell2mat(smget('K2000.V')));
end