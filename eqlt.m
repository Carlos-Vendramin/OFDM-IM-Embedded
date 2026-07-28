%-----------------------------------------%
%               Equalizaçao               %
%-----------------------------------------%
function [bits_rx, erros, x_total, y_total, y_eq] = eqlt(b, N, x_bloco, y_bloco, H_eff, Catalogo_X, Catalogo_Bits, bits_tx, erros, x_total, y_total, y_eq)

    VetorBeg = (b-1) * N + 1;
    VetorEnd = b * N;

    x_total(VetorBeg:VetorEnd) = x_bloco;
    y_total(VetorBeg:VetorEnd) = y_bloco;
    y_eq(VetorBeg:VetorEnd) = y_bloco ./ H_eff; % Apenas para visualização sem anéis

    Catalogo_X_filt = Catalogo_X .* H_eff;
    distancias = sum(abs(y_bloco - Catalogo_X_filt).^2, 1);

    %encontrando o indice de menor distancia
    [~, indice_menor] = min(distancias);
    bits_rx = Catalogo_Bits(indice_menor, :);

    erros_bloco = sum(bits_tx ~= bits_rx);
    
    erros = erros + erros_bloco;

end