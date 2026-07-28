%-------------------------------------------%
%            Transmissao                    %
%-------------------------------------------%

function [x_bloco, y_bloco, bits_tx, txstate, rxstate, tempo] = TxRx_Cascata(fs, Fator_Potencia, p, SNR_dB, LUT, N, nCicP, p1, M, K, a_filt, b_filt, txstate, rxstate, tempo, cfo_hz, Delta_tx, Delta_rx)

    bits_tx = randi([0 1], 1, p);
    bits_indice  = bits_tx(1:p1);
    bits_simbolo = bits_tx(p1+1:end); 

    dec_indice = bit2int(bits_indice.', p1);
    indices_ativos = LUT(dec_indice + 1, :); % +1 pois o MATLAB começa no índice 1

    b_sim_matriz = reshape(bits_simbolo, log2(M), K);
    simbolos_dec = bit2int(b_sim_matriz, log2(M));
    simbolos_qam = qammod(simbolos_dec, M, 'UnitAveragePower', true) * Fator_Potencia;

    x_bloco = zeros(N, 1);
    x_bloco(indices_ativos) = simbolos_qam; 

    x_tempo = ifft(x_bloco, N);
    %inserçao do prefixo ciclico
    x_Cic = [x_tempo(end-nCicP+1:end); x_tempo];

    %-------------------------------------------%
    %              Imperfeiçoes                 %
    %-------------------------------------------%
    tx_re = round(real(x_Cic)/Delta_tx)*Delta_tx;
    tx_im = round(imag(x_Cic)/Delta_tx)*Delta_tx;
    x_Cic_Quantizado = tx_re + (1i* tx_im);
    
    [x_Cic_Filtrado, txstate] = filter(b_filt, a_filt, x_Cic_Quantizado, txstate);

    % Aplicando o desvio de frequencia no canal quantizado e filtrado (CFO)
    Vetor_tempo = tempo + (0:length(x_Cic_Filtrado)-1).'/fs; 
    tempo = Vetor_tempo(end) + (1/fs); 
    
    Cfo = exp(1i*2*pi*cfo_hz*Vetor_tempo);
    x_Cic_Final = x_Cic_Filtrado .* Cfo;

    %Desse modo, o ruido atua sobre o sinal estendido pelo prefixo e quantizado
    y_Cic = awgn(x_Cic_Final, SNR_dB, 'measured');

    [rxstate, y_bloco, ~] = reconstrutor(cfo_hz, N, Vetor_tempo, y_Cic, a_filt, b_filt, Delta_rx, nCicP, rxstate);

end