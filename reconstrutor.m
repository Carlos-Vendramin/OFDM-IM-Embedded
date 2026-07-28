%--------------------------------------------%
%                  Receptor                  %
%--------------------------------------------%
function [rxstate, y_bloco, y_Cic_Final] = reconstrutor(cfo_hz, N, Vetor_tempo, y_Cic, a_filt, b_filt, Delta_rx, nCicP, rxstate)

    %compensando o CFO antes da filtragem evita rotaçao de fase por atraso de grupo
    Cfo_compensacao = exp(-1i * 2 * pi * cfo_hz * Vetor_tempo);
    y_Cic_Compensado = y_Cic .* Cfo_compensacao;

    %--------------------------------------------%
    %              Reconstruçao                  %
    %--------------------------------------------%
    [y_Cic_Filtrado, rxstate] = filter(b_filt, a_filt, y_Cic_Compensado, rxstate);

    rx_re = round(real(y_Cic_Filtrado)/Delta_rx)*Delta_rx;
    rx_im = round(imag(y_Cic_Filtrado)/Delta_rx)*Delta_rx;
    y_Cic_Final = rx_re + (1i*rx_im);

    y_tempo = y_Cic_Final(nCicP+1:end);
    y_bloco = fft(y_tempo, N);

end