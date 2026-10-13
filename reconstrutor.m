%--------------------------------------------%
%                  Receptor                  %
%--------------------------------------------%
function [rxstate, y_bloco, y_Cic_Final] = reconstrutor(cfo_hz, N, Vetor_tempo, y_Cic, a_filt, b_filt, Delta_rx, nCicP, rxstate, Compensa_CFO, Habilita_Quantizacao, Habilita_Filtro)
% Segue a ordem de um hardware real. Filtro -> Quantizaçao -> CFO
    % 1. Filtragem (Modelando o filtro anti-aliasing analógico do hardware)
    if Habilita_Filtro
        [y_Cic_Filtrado, rxstate] = filter(b_filt, a_filt, y_Cic, rxstate);
    else
        y_Cic_Filtrado = y_Cic;
    end
    
    % 2. Quantização (Modelando a amostragem pelo ADC de 12 bits)
    if Habilita_Quantizacao
        rx_re = round(real(y_Cic_Filtrado)/Delta_rx)*Delta_rx;
        rx_im = round(imag(y_Cic_Filtrado)/Delta_rx)*Delta_rx;
        y_Cic_ADC = rx_re + (1i*rx_im);
    else
        y_Cic_ADC = y_Cic_Filtrado;
    end

    % 3. Compensação de CFO Digital
    if Compensa_CFO
        Cfo_compensacao = exp(-1i * 2 * pi * cfo_hz * Vetor_tempo);
        y_Cic_Final = y_Cic_ADC .* Cfo_compensacao;
    else
        y_Cic_Final = y_Cic_ADC;
    end

    % Remoção do CP e passagem para o domínio da frequência
    y_tempo = y_Cic_Final(nCicP+1:end);
    y_bloco = fft(y_tempo, N);

end