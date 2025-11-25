function plotComparisonFor3Methods(methodStructs, methodA, methodB, methodC, outputFolder, fileSuffix)
    labelA = strrep(methodA, '_', ' ');
    labelB = strrep(methodB, '_', ' ');
    labelC = strrep(methodC, '_', ' ');
    
    if all(isfield(methodStructs, {methodA, methodB, methodC}))
        disp(['Running 3-method comparison: ' labelA ' vs ' labelB ' vs ' labelC ' ...']);

        % Extract maps
        D_A = methodStructs.(methodA).DoLP;
        D_B = methodStructs.(methodB).DoLP;
        D_C = methodStructs.(methodC).DoLP;

        A_A = methodStructs.(methodA).AoLP;
        A_B = methodStructs.(methodB).AoLP;
        A_C = methodStructs.(methodC).AoLP;

        % DoLP comparison
        compare_methods_DoLP(D_A, labelA, D_B, labelB, D_C, labelC);
        saveas(gcf, fullfile(outputFolder, sprintf('%s_%s_DoLP_%s.png', fileSuffix)));

        % AoLP comparison
        compare_methods_AoLP_hsi(A_A, labelA, A_B, labelB, A_C, labelC);
        saveas(gcf, fullfile(outputFolder, sprintf('%s_%s_AoLP_%s.png', fileSuffix)));

    else
        disp(['Skipping 3-method comparison (Missing methods: ', labelA ', ' labelB ', ' labelC ')']);
    end
end