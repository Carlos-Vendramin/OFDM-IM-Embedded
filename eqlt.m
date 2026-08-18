%-----------------------------------------%
%               Equalizaçao               %
%-----------------------------------------%
function [bits_rx, erros, x_total, y_total, y_eq] = eqlt(b, N, n, G, p_grupo, x_bloco, y_bloco, H_eff, Catalogo_X, Catalogo_Bits, bits_tx, erros, x_total, y_total, y_eq)

    VetorBeg = (b-1) * N + 1;
    VetorEnd = b * N;

    x_total(VetorBeg:VetorEnd) = x_bloco;
    y_total(VetorBeg:VetorEnd) = y_bloco;
    y_eq(VetorBeg:VetorEnd) = y_bloco ./ H_eff; % Apenas para visualização sem anéis

    bits_rx = zeros(1, p_grupo * G);
    
    for g = 1:G
        idx_start = (g-1)*n + 1;
        idx_end = g*n;
        
        y_grupo = y_bloco(idx_start:idx_end);
        H_eff_grupo = H_eff(idx_start:idx_end);
        
        Catalogo_X_filt = Catalogo_X .* H_eff_grupo;
        distancias = sum(abs(y_grupo - Catalogo_X_filt).^2, 1);

        %encontrando o indice de menor distancia
        [~, indice_menor] = min(distancias);
        bits_rx((g-1)*p_grupo + 1 : g*p_grupo) = Catalogo_Bits(indice_menor, :);
    end

    erros_bloco = sum(bits_tx ~= bits_rx);
    erros = erros + erros_bloco;
end