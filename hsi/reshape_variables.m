function [S0, S1, S2,DoLP, AoLP] = reshape_variables(S0, S1, S2, DoLP, AoLP, rows, cols) %reshape computed variables into the wanted dimensions (2d)
    S0 = reshape(S0, rows, cols);
    S1 = reshape(S1, rows, cols);
    S2 = reshape(S2, rows, cols);
    
    DoLP = reshape(DoLP, rows, cols);
    AoLP = rad2deg(reshape(AoLP, rows, cols)); % Convert to degrees

end