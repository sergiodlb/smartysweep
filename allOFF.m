function allOFF()
    fprintf('all OFF..!\n');
%     list = {'Vdd', 'Vg', 'Vbg', 'Vdc'};
    list = {'Vdd', 'Vg', 'Vbg', 'Vdc', 'ZI.OUT1AMP1', 'ZI.OUT2AMP2'};
    smset(list, 0, 1);
    fprintf('%s\n', sprintf('%s\t', list{:}));
    fprintf('%s\n', sprintf('%.4g\t', cell2mat(smget(list))));
end