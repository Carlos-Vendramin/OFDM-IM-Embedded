clear all; clc; close all;
%% Configuraçao do transmissor OFDM-IM
N = 16;              % total de subportadoras
K = 2;              % Número de subportadoras ativas
M = 16;              % Ordem da modulação (M-QAM)
num_blocos = 10000;  % Numero de amostras a serem transmitidas no loop
%nFFT = 124; %Numero de FFT'S
nCicP = 9;


EbNo_vec = 0:2:20;
BER = zeros(1, length(EbNo_vec));

%----------------------------------%
%       Correçao da potencia       %
%----------------------------------%

Pot_N=N/K;
Fator_Potencia=sqrt(Pot_N);


%-------------------------------------------%
%       Frequencia e largura de banda       %
%-------------------------------------------%
[Fc, Bw, fs, F_Corte, a_filt,b_filt, H_freq, w, H_eff] =freq_bw(N,2.4e9, 20e6, 20e6, 9.9e6);
%Correto
%----------------------------------


%-------------------------------------------%
%              Modelagem CFO                %
%-------------------------------------------%
%Carrier frequency offset
%[cfo_hz] = cfo(Crystal_ppm, Fc)
[cfo_hz] = cfo(10e-6, Fc);
%----------------------------------



%----------------------------------------------------%
%            Declaraçao de vetores RX/TX             %
%----------------------------------------------------%
[x_total, y_total, y_eq] = rxtx_Vec(N, num_blocos);
%correto
%----------------------------------

%-------------------------------------------%
%            Lookup table (LUT)             %
%-------------------------------------------%
[p1, p2, p, combinacoes_LUT, LUT] = LUT(N,K,M);
%-----------------------------------------------------------------------------


%-----------------------------------------------------------------------------%
%           Catalogo de verossimilhança (Maximum Likelihood - ML)             %
%-----------------------------------------------------------------------------%
combinacoes=2^p;
Catalogo_X = (zeros(N, combinacoes));
Catalogo_Bits = zeros(combinacoes, p);
for c=1:combinacoes
[Catalogo_X, qam_sum, Catalogo_Bits, b_indice, b_sim] = ML(p, c, Catalogo_X, Catalogo_Bits, LUT, p1,M,K,Fator_Potencia,N);
end
%-----------------------------------------------------------------------------


%-------------------------------------------%
%            Transmissao                    %
%-------------------------------------------%
disp('--- TRANSMISSOR ---');
disp('Fazendo um loop de transmissao, pra poder gerar uma grande amostra')
disp('Logica monte carlo para a curva de BER');

max_tx_hardware = 4.0;
max_rx_hardware = 4.0;
bits_adc = 12;
bits_dac = 8;
Delta_tx = 2 * max_tx_hardware / (2^bits_dac);
Delta_rx = 2 * max_rx_hardware / (2^bits_adc);
 
for snr_indice = 1:length(EbNo_vec)
    
    EbNo_dB = EbNo_vec(snr_indice);
    SNR_dB = EbNo_dB + 10*log10(p/(N+nCicP));
    
    erros = 0; 
    tempo = 0;
    txstate = zeros(max(length(a_filt)-1, length(b_filt)-1), 1);
    rxstate = zeros(max(length(a_filt)-1, length(b_filt)-1), 1);
    

%--------------------------------------------------------%
%               variaveis sensoriamento                  %
%--------------------------------------------------------%

distancia=1;
velocidade=1;

Radar_tx_Matriz = zeros(N, num_blocos);
Radar_rx_Matriz = zeros(N, num_blocos);




%----------------------------------------------%

    for b = 1:num_blocos
        
        [x_bloco, y_bloco, bits_tx, txstate, rxstate, tempo] = TxRx_Cascata(fs, Fator_Potencia, p, SNR_dB, LUT, N, nCicP, p1, M, K, a_filt, b_filt, txstate, rxstate, tempo, cfo_hz, Delta_tx, Delta_rx);
        
        [bits_rx, erros, x_total, y_total, y_eq] = eqlt(b, N, x_bloco, y_bloco, H_eff, Catalogo_X, Catalogo_Bits,bits_tx, erros, x_total, y_total, y_eq);
            
    

    %-----------------------------------------%
    %              Sensoriamento              %
    %-----------------------------------------%
    
        fase_range=exp(-1i * 2 * pi * (0:N-1).' * distancia / N);
        fase_doppler= exp(1i * 2 * pi * velocidade * b / num_blocos);
    
        y_radar_bloco = (x_bloco .* fase_range) * fase_doppler;
        y_radar_bloco = awgn(y_radar_bloco, SNR_dB, 'measured');
        
        Radar_tx_Matriz(:, b) = x_bloco;
        Radar_rx_Matriz(:, b) = y_radar_bloco;
    %-----------------------------------------%

    end
    
    filtro_casado = Radar_rx_Matriz .* conj(Radar_tx_Matriz);
    Range = ifft(filtro_casado, N, 1);
    Doppler = fftshift(fft(Range, num_blocos, 2), 2);


    %-----------------------------------------%
    %      Iniciando a calculo da BER         %
    %-----------------------------------------%
    total_bits_TX = num_blocos * p;
    BER(snr_indice) = erros / total_bits_TX;
    
    disp(['SNR = ', num2str(SNR_dB), ' dB | BER = ', num2str(BER(snr_indice))]);
end


%-----------------------------------------%
%               Plotagem                  %
%-----------------------------------------%

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));

fig_constelacao=figure;
%pontos recebidos com ruido (vermelho)
plot(real(y_eq), imag(y_eq), 'r.', 'MarkerSize', 8);
hold on;
%transmitidos e ideiais (circulos pretos)
plot(real(x_total), imag(x_total), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

%formataçao do grafico
grid on;
xlabel('Fase (I) - Parte Real');
ylabel('Quadratura (Q) - Parte Imaginária');
title(['Constelação OFDM-IM (N=', num2str(N), ', K=', num2str(K), ', ', num2str(M), '-QAM) - SNR = ', num2str(SNR_dB), ' dB']);
legend('Sinal Recebido (com ruído)', 'Sinal Transmitido (Ideal)', 'Location', 'best');

% Fixa os limites dos eixos X e Y para -2 a +2 para visualização simétrica
axis([-4 4 -4 4]);
disp('--- RECEPTOR ---');
disp('Sinal recebido com ruído na Frequência:');
disp(y_bloco);
png_constelacao = ['constelacao_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_constelacao, png_constelacao, 'Resolution', 300);


% 4. Curva BER
fig_BER=figure;
% semilogy é usado porque a BER decai exponencialmente (eixo Y logarítmico)
semilogy(EbNo_vec, BER, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(['Desempenho de BER do OFDM-IM (N=', num2str(N), ', K=', num2str(K), ', ', num2str(M), '-QAM)'])
axis([min(EbNo_vec) max(EbNo_vec) 1e-5 1]); % Fixa os limites do gráfico
png_BER = ['BER_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_BER, png_BER, 'Resolution', 300);

%-----------------------------------------------------------------------------

% 5. Mapa range doppler
fig_RangeDoppler=figure;
doppler_axis = linspace(0, 10, num_blocos);
range_axis = 0:(N-1);

mesh(doppler_axis, range_axis, 10*log10(abs(Range).^2));
view(2);
xlabel('Frequência Doppler (Velocidade)');
ylabel('Range (Distância)');
title(['Mapa Range-Doppler OFDM-IM | SNR = ', num2str(SNR_dB), ' dB']);
colorbar;
png_RD = ['RangeDoppler_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_RangeDoppler, png_RD, 'Resolution', 300);