%-------------------------------------------%
%            Transmissao                    %
%-------------------------------------------%

function [x_bloco, y_bloco, bits_tx, txstate, rxstate, tempo] = TxRx_Cascata(fs, Fator_Potencia, p_grupo, p_total, SNR_dB, LUT, N, n, G, nCicP, p1, M, k, a_filt, b_filt, txstate, rxstate, tempo, cfo_hz, Delta_tx, Delta_rx, Habilita_CFO, Compensa_CFO, Habilita_Quantizacao, Habilita_Filtro)

    bits_tx = randi([0 1], 1, p_total);
    x_bloco = zeros(N, 1);
    
    for g = 1:G
        offset_bit = (g-1) * p_grupo;
        bits_grupo = bits_tx(offset_bit+1 : offset_bit+p_grupo);
        
        bits_indice  = bits_grupo(1:p1);
        bits_simbolo = bits_grupo(p1+1:end); 

        dec_indice = bit2int(bits_indice.', p1);
        indices_ativos = LUT(dec_indice + 1, :); % +1 pois o MATLAB começa no índice 1

        b_sim_matriz = reshape(bits_simbolo, log2(M), k);
        simbolos_dec = bit2int(b_sim_matriz, log2(M));
        simbolos_qam = qammod(simbolos_dec, M, 'UnitAveragePower', true) * Fator_Potencia;

        offset_idx = (g-1) * n;
        x_bloco(indices_ativos + offset_idx) = simbolos_qam; 
    end

    x_tempo = ifft(x_bloco, N);
    %inserçao do prefixo ciclico
    x_Cic = [x_tempo(end-nCicP+1:end); x_tempo];

    %-------------------------------------------%
    %              Imperfeiçoes                 %
    %-------------------------------------------%
    if Habilita_Quantizacao
        max_tx = (2^11) * Delta_tx; 
        re_limitado = max(min(real(x_Cic), max_tx), -max_tx);
        im_limitado = max(min(imag(x_Cic), max_tx), -max_tx);
        
        tx_re = round(re_limitado/Delta_tx)*Delta_tx;
        tx_im = round(im_limitado/Delta_tx)*Delta_tx;
        x_Cic_Quantizado = tx_re + (1i* tx_im);
    else
        x_Cic_Quantizado = x_Cic;
    end

    if Habilita_Filtro
        [x_Cic_Filtrado, txstate] = filter(b_filt, a_filt, x_Cic_Quantizado, txstate);
    else
        x_Cic_Filtrado = x_Cic_Quantizado;
    end
    
    % Aplicando o desvio de frequencia no canal quantizado e filtrado (CFO)
    Vetor_tempo = tempo + (0:length(x_Cic_Filtrado)-1).'/fs; 
    tempo = Vetor_tempo(end) + (1/fs);

   if Habilita_CFO
        Cfo = exp(1i*2*pi*cfo_hz*Vetor_tempo);
        x_Cic_Final = x_Cic_Filtrado .* Cfo;
    else
        x_Cic_Final = x_Cic_Filtrado;
    end
    
    %Desse modo, o ruido atua sobre o sinal estendido pelo prefixo e quantizado
    y_Cic = awgn(x_Cic_Final, SNR_dB, 'measured');

    % Atualize a chamada do reconstrutor passando as duas flags de CFO separadas
    [rxstate, y_bloco, ~] = reconstrutor(cfo_hz, N, Vetor_tempo, y_Cic, a_filt, b_filt, Delta_rx, nCicP, rxstate, Compensa_CFO, Habilita_Quantizacao, Habilita_Filtro);
end